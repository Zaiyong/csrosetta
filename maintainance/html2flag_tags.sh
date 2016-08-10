cat html2 | awk -v FS='</td><td>' '{print $1}' | awk -v FS='<tr><td>' '{print $2}' | awk 'NF>0 {printf("%s, ",$1)}' 

