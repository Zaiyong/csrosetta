for i in $( cat different_files | grep -v message_listening | grep -v -f excludes | sed s@vanilla/@@ ); do  cp rosetta/$i patched/$i; done
