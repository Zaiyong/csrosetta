. $HOME/csrosetta3/com/init
path=`pwd`
rm -Rf $HOME/csrosetta3/help_output
mkdir -p $HOME/csrosetta3/help_output
cd $HOME/csrosetta3/help_output

for i in $( ls ../com/* | grep -v median | grep -v RELEASE | grep -v IGNORE | grep -v -f ../com/_IGNORE | grep -f ../com/_RELEASE  ); do
    echo get raw help for $i
    $i -h > help_$( basename $i);
done

for method in $( ../com/setup_target.py -h | awk 'BEGIN {fine=1} /examples:/ {fine=0} fine==1 {print}' | grep rasrec | grep -v 'setup_target.py' | grep -v '\[' | sed s/,/\ /g | sed s/-method// | sed 's/{//' | sed 's/}//' | tr -cs "[:alpha:]" "\n" | grep -f ../flag_library/methods/_RELEASE ); do
    echo get method specific options for $method...
    ../com/setup_target.py -h -method $method > method_spec_help_setup_target_"$method".py
#filter out only the method-specific part...
#    cat method_spec_help_setup_target_"$method".py | awk 'BEGIN {no=2} $1=="method-options:" {no=1;next} NF>1 && no<=0 {print; next} no<2 && NF>1 {no=no-1; next}' > method_spec_help_setup_target_"$method"_extracted.py
    cat method_spec_help_setup_target_"$method".py | $path/argparse2html.py -method $method > help_setup_target_$method.html
done


for method in $( ../com/setup_run.py -h | awk 'BEGIN {fine=1} /examples:/ {fine=0} fine==1 {print}' | grep rasrec | grep -v 'setup_target.py' | grep -v '\[' | sed s/,/\ /g | sed s/-method// | sed 's/{//' | sed 's/}//' | tr -cs "[:alpha:]" "\n" | grep -f ../flag_library/methods/_RELEASE ); do
    echo get method specific options for $method...
    ../com/setup_run.py -h -method $method > method_spec_help_setup_run_"$method".py
#filter out only the method-specific part...
#    cat method_spec_help_setup_target_"$method".py | awk 'BEGIN {no=2} $1=="method-options:" {no=1;next} NF>1 && no<=0 {print; next} no<2 && NF>1 {no=no-1; next}' > method_spec_help_setup_target_"$method"_extracted.py
    cat method_spec_help_setup_run_"$method".py | $path/argparse2html.py -method $method > help_setup_run_$method.html
done

find . -size 0 -exec rm {} \;
grep usage help* | awk -v FS=":" '{print $1}' > proper_help_files.txt
for i in $( cat proper_help_files.txt | grep -v ~ ); do cat $i | $path/argparse2html.py > $( basename $i .py ).html; done
rm *~*html
cat <<EOF > upload_application_help.xml
<?xml version="1.0" encoding="UTF-8" ?>
<rss version="2.0">
<channel>
        <title>RSS Title</title>
EOF
ff=upload_application_help.xml
for i in $( ls *html ); do
    APP=$( echo $i | sed s/help_// | sed s/.html// )
    #take away the setup_target_ from specific method option help
    TIT=$( echo $APP | sed 's/setup_target_/choice: /' | sed 's/setup_run_/choice: /' )
    echo '<item>' >> $ff
    echo '<title>'$TIT'</title>' >> $ff
    echo '<guid>'$APP'</guid>' >> $ff
    echo '<description>' >> $ff
    cat $i | awk '{printf("<![CDATA[%s]]>\n",$0)}' >> $ff
    echo '</description>' >> $ff
    echo '</item>' >> $ff
done
cat <<EOF >> $ff
</channel>
</rss>
EOF

scp upload_application_help.xml di34ron@webdev1.lrz.de:/nfs/web_tum/www/n/di34ron/webserver/drupal/sites/default/files/feeds/upload_application_help_0.xml
echo 'data has been uploaded to webserver'
echo 'now go to www.csrosetta.org/import/import_application_pages'
cd -
