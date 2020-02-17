#!/bin/csh

echo "Yo, move the data to ? (y/n): "
set order = $<
if ($order == "y") then
   mv gal_*.parts /scratch/shuiyao/scidata/galparts/p50n288zw
else
   echo "Abort then. No big deal."
endif

