use 5.010000;
use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.

WriteMakefile(
    NAME              => 'MOODS',
    VERSION_FROM      => 'lib/MOODS.pm', # finds $VERSION
    PREREQ_PM         => {}, # e.g., Module::Name => 1.1
    ($] >= 5.005 ?     ## Add these new keywords supported since 5.005
      (ABSTRACT_FROM  => 'lib/MOODS.pm', # retrieve abstract from module
       AUTHOR         => 'Petri J Martinmaki <pmartinm@>') : ()),
    LIBS              => ['-L../src -lpssm '], # e.g., '-lm'
    DEFINE            => '', # e.g., '-DHAVE_SOMETHING'
    INC               => '-I. -I../src', # e.g., '-I. -I/usr/include/other'
    CC                => "g++",
    LD                => "g++",  
	# Un-comment this if you add C files to link with later:
    #OBJECT            => '$(O_FILES)', # link all the C files too
);
