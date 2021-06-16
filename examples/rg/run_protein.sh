#!/bin/bash
python rg.py ../data/SARS_CoV_2_protease.gro ../data/SARS_CoV_2_protease.xtc protein \
--ocsv out/protease_all.csv --oimage out/protease_all.png

python rg.py ../data/SARS_CoV_2_protease.gro ../data/SARS_CoV_2_protease.xtc heavy \
--ocsv out/protease_heavy.csv --oimage out/protease_heavy.png

python rg.py ../data/SARS_CoV_2_protease.gro ../data/SARS_CoV_2_protease.xtc backbone \
--ocsv out/protease_backbone.csv --oimage out/protease_backbone.png

python rg.py ../data/SARS_CoV_2_protease.gro ../data/SARS_CoV_2_protease.xtc alpha_C \
--ocsv out/protease_alpha_C.csv --oimage out/protease_alpha_C.png
