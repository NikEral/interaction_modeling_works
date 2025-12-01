#!/bin/bash

for i in {1..5}
do 
echo "$i"
	root -l -q "script.C(10, $i*(-5), $i*5, $i*10, \"example$i.root\")"
done
