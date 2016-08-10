mkdir -p $csrosettaDir/com_alias; cd $csrosettaDir/com_alias; for i in $( ls ../com/*py ../com/*com ../com/*.pl ); do ln -sf $i $( basename $i .py ); done


