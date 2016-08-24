############################################
#
#   Beat It     0.1
#
############################################

Install instructions:

Setup your project folder:

    BeatIt/
        beatit  (git repository)
        build   (build folder)
        install (prefix folder)

Copy the configure.sh file  from the git repository in your build folder.

BeatIt/build$ cp ../beatit/do-configure.sh my-configure.sh

Edit the my-configure.sh file to set your own libraries

From your build folder run
BeatIt/build$ sh my-configure.sh
