#!/bin/bash
species=$1
filename=$2
grep -C 10 "<Hit_num>1</Hit_num>" $filename | grep -B 10 $species | grep "<Iteration_query-def>"
