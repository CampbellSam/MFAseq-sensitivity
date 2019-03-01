#!/bin/bash


for bed in *.bed; do

echo "no window $bed"
python ~/git/originssimulation/readBed.py --bed $bed --chr LmjF.36 --start 155190 --end 156279
python ~/git/originssimulation/readBed.py --bed $bed --chr LmjF.36 --start 1110127 --end 1116528

echo "window $bed"
python ~/git/originssimulation/readBed.py --bed $bed --chr LmjF.36 --start 155190 --end 156279 --window 30000
python ~/git/originssimulation/readBed.py --bed $bed --chr LmjF.36 --start 1110127 --end 1116528 --window 30000

echo "chr $bed"
python ~/git/originssimulation/readBed.py --bed $bed --chr LmjF.36
python ~/git/originssimulation/readBed.py --bed $bed --chr LmjF.36
done