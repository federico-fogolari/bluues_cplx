rm -f exp* std* test_* trian*
bluues_cplx  pdb1ay7.pqr test_ns A B -ns test_ns_aux -srf -srfpot -pa 0.1
bluues_cplx  pdb1ay7.pqr test_msms A B -ns test_msms_aux -pa 0.1 -lite
rm -f exp* std* trian*

