
cat t | sed 's/NEW_OPT( /<tr><td>-/' | sed 's/\"/\ /' | sed 's/,/)/' | sed 's@)@</td><td>@' | sed 's/<td>\ \ /<td>/' | awk -v FS=',' '{print $1}' | sed 's@\"@</td></tr>@'

