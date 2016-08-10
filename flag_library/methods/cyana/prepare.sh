#!/bin/bash

### usually empty
mv CALC.cya CALC.cya.tmp
rand=$(sleep 10;date -- "+%s")
cat CALC.cya.tmp | sed s/210003/$rand/g > CALC.cya