svn co https://devel-ik.fzk.de/svn/mc/crmc/trunk

mv trunk/src/*.inc*         src/
mv trunk/src/*.f            src/
mv trunk/src/*.F            src/
mv trunk/src/epos/*.f       src/
mv trunk/src/sibyll/*.f     src/
mv trunk/src/qgsjet/*.f     src/
mv trunk/src/qgsjetII/*.f   src/
mv trunk/src/crmc.param.in  crmc.param
mv trunk/tabs/              .

rm src/epos-app.f 

rm -rf trunk

sed -i 's%@CRMCROOT@%${PWD}%g' crmc.param #put in the tabs location
sed -i 's%Length  1.%Length  0.02997924%g' crmc.param #put it to decay 10^-12s
sed -i '403,441 s%^%c%' src/crmc-aaa.f #comment out random generator
sed -i 's%nmxhep=9990%nmxhep=99900%g' src/crmc-aaa.f src/epos.inc #increase number of allowed particles
sed -i 's%hepevt%hepcom%g' src/epos.inc #this is a fix for a crash
