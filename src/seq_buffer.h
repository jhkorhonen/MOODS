// seq_buffer - tools for itarating over large char arrays
// Copyright (C) 2009  Petri Martinmäki
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version, or under the terms of the Biopython
// License.

#ifndef _SEQ_BUFFER_H_
#define _SEQ_BUFFER_H_

#include <cstdio>
#include <cstdlib>
#include <cassert>

//typedef unsigned long int size_type;

/**
 * Interface, whose implementation is needed by SeqIterator class.
 */
class SeqSourceI {
public:
	/**
	 * Read data must write positive char values to memory area starting from dst.
	 *
	 * @return must be a number of written characters.
	 */
	virtual int read_data(char *dst, int length) = 0;

	/**
	 * Resets the source in case one wants to start iterating over.
	 */
	virtual void reset() = 0;

	/**
	 * True if there is more data available, false otherwise.
	 */
	virtual bool eof() = 0;
};

/**
 * SeqIterator allows user to access sequence of data. SeqIterator polls and buffers
 * data from a child class of SeqSourceI.
 *
 * Usage:
 * 		SeqIterator it(source, 1024);
 * 		while(it.hasData()) {
 * 			do_something(*it);
 * 			++it;
 * 		}
 * 		it.reset();
 * 		if(it.hasData(40)) {
 * 			do_something(it[39]);
 * 		}
 */
class SeqIterator {
	SeqSourceI *seq_source;
	char *buffer;
	char *c_pos;
	char *data_end;
	char *buffer_end;
	int buffer_size;

	unsigned long moffset;
	bool renewed;

public:
	SeqIterator(SeqSourceI *source, int buffer_size) : seq_source(source), buffer_size(buffer_size) {
		buffer = (char *) malloc(buffer_size * sizeof(char));
		buffer_end = buffer + buffer_size;
		data_end = NULL;
		c_pos = NULL;
		renewed = true;
		reset();
	}
	SeqIterator(const SeqIterator &other) {
		buffer = other.buffer;
		buffer_end = other.buffer_end;
		data_end = other.data_end;
		c_pos = other.c_pos;
		renewed = other.renewed;
		seq_source = other.seq_source;
		buffer_size = other.buffer_size;
	}
	~SeqIterator() {
		free(buffer);
	}
	SeqIterator& operator ++() {
		c_pos++;
		return *this;
	}
	SeqIterator& operator +=(int value) {
		c_pos += value;
		return *this;
	}
	char &operator *() {
		//assert(c_pos >= buffer);
		//assert(c_pos < buffer_end);
		return *c_pos;
	}
	char &operator [](int pos) {
		//assert(c_pos + pos >= buffer);
		//assert(c_pos + pos < buffer_end);
		return *(c_pos +pos);
	}

	bool hasData() {
		return c_pos < data_end ? true : renew();
	}
	bool hasData(int length) {
		if(c_pos + length <= data_end) {
			return true;
		}
		else {
			while(renew() && c_pos + length <= data_end);
			return c_pos + length <= data_end;
		}
	}
	int buffered() {
		return data_end - c_pos;
	}
	bool buffered(int m) {
		return data_end - c_pos >= m;
	}
	void reset() {
		c_pos = buffer;
		moffset = 0;
		if(renewed) {
			seq_source->reset();
			data_end = buffer + seq_source->read_data(buffer, buffer_end - buffer);
			renewed = false;
		}
	}
	unsigned long position() {
		return moffset + c_pos - buffer;
	}
	bool readOffset(int length) {
		//printf("Offset read asked %d\n", length);
		int readsum = 0;
		int readlen;
		int tmp = buffer_size;
		while(readsum < length) {
			if(length - readsum < buffer_size) {
				tmp = length - readsum;
			}
			else
				tmp = buffer_size;
			readlen = seq_source->read_data(buffer, tmp);
			if(readlen <= 0)
				return false;
			readsum += readlen;
			//printf("Offset read: %d\n", readlen);
		}
		return true;
	}
	bool renew() {
		if(seq_source->eof() || (c_pos == buffer && data_end == buffer_end)) {
			return false;
		}
		moffset += c_pos - buffer;
		renewed = true;
		if(c_pos == data_end) {
			c_pos = buffer;
			data_end = buffer;
		}
		else if(c_pos > data_end) {
			if(!readOffset(c_pos - data_end)) {
				return false;
			}
			c_pos = buffer;
			data_end = buffer;
		}
		else {
			char *tmp = buffer;
			while(c_pos != data_end) {
				*tmp = *c_pos;
				c_pos ++;
				tmp ++;
			}
			c_pos = buffer;
			data_end = tmp;
		}
		data_end = data_end + seq_source->read_data(data_end, buffer_end - data_end);
		if(c_pos < data_end) {
			return true;
		}
		else
			return false;
	}
};

#endif
