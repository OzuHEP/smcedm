#!/bin/bash
#cd /mnt/harddisk2/wjets/wellnu01234j_5f_LO_MLM/Events/run_04_$1/
#pigz -d tag_1_pythia8_events.hepmc.gz

cd /home/bora/softwares/MG5_aMC_v2_6_7/Delphes
./DelphesHepMC /usr/local/share/delphes_cards/delphes_card_CMS_skimmed.tcl /mnt/harddisk2/wjets/wellnu01234j_5f_LO_MLM/Events/run_04_$1/delphes_$1.root /mnt/harddisk2/wjets/wellnu01234j_5f_LO_MLM/Events/run_04_$1/wjets_pythia8_events.hepmc

cd /mnt/harddisk2/wjets/wellnu01234j_5f_LO_MLM/Events/run_02_$1/
pigz -1 wjets_pythia8_events.hepmc

