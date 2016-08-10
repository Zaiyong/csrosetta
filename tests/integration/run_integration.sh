#!/bin/bash

source ../../com/init
rm -r new
commit=$(head ../../.git/RENAMED-REF -n 1 | awk '{print $1}' )
ex_commit=$(head ../../.git/ORIG_HEAD -n 1 | awk '{print $1}' )
author=$( git log -1 | grep Author | awk '{print $2}' )
date=$(git log -1 | grep Date | awk '{print $2,$3,$4,$5,$6}')
python2.7 ./integration.py --log=integration_"$commit".log --cur_commit=$commit --ex_commit=$ex_commit
mkdir -p integration_log
mv integration_"$commit".log integration_log/
rm ref

if [ ! -e ref_"$commit" ];then
mv new ref_"$commit"
fi
ln -s ref_"$commit" ref

mv buildbot_show.html buildbot_old.html
cat buildbot_old.html | head -n "$NF"-3 > buildbot_show.html
echo "<tr>" >> buildbot_show.html
echo "<td>$commit</td>" >> buildbot_show.html
echo "<td>$author</td>" >> buildbot_show.html
echo "<td>$date</td>" >> buildbot_show.html
fail=$( cat integration_log/integration_"$commit".log | grep FAIL | wc -l)

if [ $fail -gt 0 ];then
	bg_color=red
	r_txt=failure
else
	bg_color=green
	r_txt=success
fi
cur_dir=$(pwd)
link=integration_log/integration_"$commit".log

echo "<td style=\"background-color:$bg_color;\">" >> buildbot_show.html
echo "<a href=\"$link\">" >> buildbot_show.html
echo "$r_txt </a></td>" >> buildbot_show.html

cd ..
mkdir -p unittest_log
cp unittest.log unittest_log/unittest_"$commit".log
chmod -R a+r+x unittest_log

pass=$(cat unittest.log | grep FAILED | wc -l )


if [ $pass -gt 0 ];then
	u_bg_color=red
	u_r_txt=failure
else
	u_bg_color=green
	u_r_txt=success
fi
cd -
u_link=../unittest_log/unittest_"$commit".log

echo "<td style=\"background-color:$u_bg_color;\">" >> buildbot_show.html
echo "<a href=\"$u_link\">" >> buildbot_show.html
echo "$u_r_txt </a></td>" >> buildbot_show.html

echo "</tr>" >> buildbot_show.html
echo "</table>" >> buildbot_show.html
echo "</body>" >> buildbot_show.html
echo "</html>" >> buildbot_show.html

chmod -R a+r+x *
