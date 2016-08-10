#!/bin/sh

site=http://www.csr-dev.org/
name=shell_user
pass=csrosetta

cookies=/tmp/cookie.txt
wget -O /dev/null --save-cookies /tmp/ba-cookies.txt --keep-session-cookies --load-cookies $cookies "${site}user"
wget --keep-session-cookies --save-cookies $cookies --load-cookies $cookies -O /dev/null \
        --post-data="name=$name&pass=$pass&op=Log%20in&form_id=user_login" \
        "${site}user?destination=login_redirect"
rm -f vall toolbox
wget --keep-session-cookies --save-cookies $cookies --load-cookies $cookies "${site}/downloads/vall"
link=$( grep '\<tr><td><a' vall | sed 's/href=/\ /' | tr '"' ' ' | awk '{print $2}')
wget -nc --keep-session-cookies --save-cookies $cookies --load-cookies $cookies "${site}/$link"

wget --keep-session-cookies --save-cookies $cookies --load-cookies $cookies "${site}/downloads/toolbox"
link=$( grep '\<tr><td><a' toolbox | sed 's/href=/\ /' | tr '"' ' ' | awk '{print $2}')
wget -nc --keep-session-cookies --save-cookies $cookies --load-cookies $cookies "${site}/$link"
