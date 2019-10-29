\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  0 z0 + z1 + 0 z2 + 0 z3 + z4 + z5 + 0 z6 + 0 z7 + 0 z8 + 0 z9
Subject To
 root_flow: - f081 + f801 - f381 + f831 - f581 + f851 - f681 + f861 + f891
   - f981 = 1
 root_flow: - f084 + f804 - f384 + f834 - f584 + f854 - f684 + f864 + f894
   - f984 = 1
 root_flow: - f085 + f805 - f385 + f835 - f585 + f855 - f685 + f865 + f895
   - f985 = 1
 root_flow: - f080 + f800 - f380 + f830 - f580 + f850 - f680 + f860 + f890
   - f980 = 1
 root_flow: - f082 + f802 - f382 + f832 - f582 + f852 - f682 + f862 + f892
   - f982 = 1
 root_flow: - f083 + f803 - f383 + f833 - f583 + f853 - f683 + f863 + f893
   - f983 = 1
 root_flow: - f086 + f806 - f386 + f836 - f586 + f856 - f686 + f866 + f896
   - f986 = 1
 root_flow: - f087 + f807 - f387 + f837 - f587 + f857 - f687 + f867 + f897
   - f987 = 1
 root_flow: - f089 + f809 - f389 + f839 - f589 + f859 - f689 + f869 + f899
   - f989 = 1
 flow_conservation: - f031 + f301 - f051 + f501 - f061 + f601 - f081 + f801
   = 0
 flow_conservation: = 0
 flow_conservation: f031 - f301 - f341 + f431 - f351 + f531 - f371 + f731
   - f381 + f831 - f391 + f931 = 0
 flow_conservation: f141 - f411 + f341 - f431 - f451 + f541 - f471 + f741
   - f491 + f941 = 0
 flow_conservation: f051 - f501 + f151 - f511 + f351 - f531 + f451 - f541
   - f561 + f651 - f571 + f751 - f581 + f851 - f591 + f951 = 0
 flow_conservation: f061 - f601 + f561 - f651 - f681 + f861 = 0
 flow_conservation: f171 - f711 + f371 - f731 + f471 - f741 + f571 - f751
   - f791 + f971 = 0
 flow_conservation: f191 - f911 + f391 - f931 + f491 - f941 + f591 - f951
   + f791 - f971 + f891 - f981 = 0
 flow_conservation: - f034 + f304 - f054 + f504 - f064 + f604 - f084 + f804
   = 0
 flow_conservation: - f144 + f414 - f154 + f514 - f174 + f714 - f194 + f914
   = 0
 flow_conservation: = 0
 flow_conservation: f034 - f304 - f344 + f434 - f354 + f534 - f374 + f734
   - f384 + f834 - f394 + f934 = 0
 flow_conservation: f054 - f504 + f154 - f514 + f354 - f534 + f454 - f544
   - f564 + f654 - f574 + f754 - f584 + f854 - f594 + f954 = 0
 flow_conservation: f064 - f604 + f564 - f654 - f684 + f864 = 0
 flow_conservation: f174 - f714 + f374 - f734 + f474 - f744 + f574 - f754
   - f794 + f974 = 0
 flow_conservation: f194 - f914 + f394 - f934 + f494 - f944 + f594 - f954
   + f794 - f974 + f894 - f984 = 0
 flow_conservation: - f035 + f305 - f055 + f505 - f065 + f605 - f085 + f805
   = 0
 flow_conservation: - f145 + f415 - f155 + f515 - f175 + f715 - f195 + f915
   = 0
 flow_conservation: = 0
 flow_conservation: f035 - f305 - f345 + f435 - f355 + f535 - f375 + f735
   - f385 + f835 - f395 + f935 = 0
 flow_conservation: f145 - f415 + f345 - f435 - f455 + f545 - f475 + f745
   - f495 + f945 = 0
 flow_conservation: f065 - f605 + f565 - f655 - f685 + f865 = 0
 flow_conservation: f175 - f715 + f375 - f735 + f475 - f745 + f575 - f755
   - f795 + f975 = 0
 flow_conservation: f195 - f915 + f395 - f935 + f495 - f945 + f595 - f955
   + f795 - f975 + f895 - f985 = 0
 flow_conservation: - f140 + f410 - f150 + f510 - f170 + f710 - f190 + f910
   = 0
 flow_conservation: = 0
 flow_conservation: f030 - f300 - f340 + f430 - f350 + f530 - f370 + f730
   - f380 + f830 - f390 + f930 = 0
 flow_conservation: f140 - f410 + f340 - f430 - f450 + f540 - f470 + f740
   - f490 + f940 = 0
 flow_conservation: f050 - f500 + f150 - f510 + f350 - f530 + f450 - f540
   - f560 + f650 - f570 + f750 - f580 + f850 - f590 + f950 = 0
 flow_conservation: f060 - f600 + f560 - f650 - f680 + f860 = 0
 flow_conservation: f170 - f710 + f370 - f730 + f470 - f740 + f570 - f750
   - f790 + f970 = 0
 flow_conservation: f190 - f910 + f390 - f930 + f490 - f940 + f590 - f950
   + f790 - f970 + f890 - f980 = 0
 flow_conservation: - f032 + f302 - f052 + f502 - f062 + f602 - f082 + f802
   = 0
 flow_conservation: - f142 + f412 - f152 + f512 - f172 + f712 - f192 + f912
   = 0
 flow_conservation: f032 - f302 - f342 + f432 - f352 + f532 - f372 + f732
   - f382 + f832 - f392 + f932 = 0
 flow_conservation: f142 - f412 + f342 - f432 - f452 + f542 - f472 + f742
   - f492 + f942 = 0
 flow_conservation: f052 - f502 + f152 - f512 + f352 - f532 + f452 - f542
   - f562 + f652 - f572 + f752 - f582 + f852 - f592 + f952 = 0
 flow_conservation: f062 - f602 + f562 - f652 - f682 + f862 = 0
 flow_conservation: f172 - f712 + f372 - f732 + f472 - f742 + f572 - f752
   - f792 + f972 = 0
 flow_conservation: f192 - f912 + f392 - f932 + f492 - f942 + f592 - f952
   + f792 - f972 + f892 - f982 = 0
 flow_conservation: - f033 + f303 - f053 + f503 - f063 + f603 - f083 + f803
   = 0
 flow_conservation: - f143 + f413 - f153 + f513 - f173 + f713 - f193 + f913
   = 0
 flow_conservation: = 0
 flow_conservation: f143 - f413 + f343 - f433 - f453 + f543 - f473 + f743
   - f493 + f943 = 0
 flow_conservation: f053 - f503 + f153 - f513 + f353 - f533 + f453 - f543
   - f563 + f653 - f573 + f753 - f583 + f853 - f593 + f953 = 0
 flow_conservation: f063 - f603 + f563 - f653 - f683 + f863 = 0
 flow_conservation: f173 - f713 + f373 - f733 + f473 - f743 + f573 - f753
   - f793 + f973 = 0
 flow_conservation: f193 - f913 + f393 - f933 + f493 - f943 + f593 - f953
   + f793 - f973 + f893 - f983 = 0
 flow_conservation: - f036 + f306 - f056 + f506 - f066 + f606 - f086 + f806
   = 0
 flow_conservation: - f146 + f416 - f156 + f516 - f176 + f716 - f196 + f916
   = 0
 flow_conservation: = 0
 flow_conservation: f036 - f306 - f346 + f436 - f356 + f536 - f376 + f736
   - f386 + f836 - f396 + f936 = 0
 flow_conservation: f146 - f416 + f346 - f436 - f456 + f546 - f476 + f746
   - f496 + f946 = 0
 flow_conservation: f056 - f506 + f156 - f516 + f356 - f536 + f456 - f546
   - f566 + f656 - f576 + f756 - f586 + f856 - f596 + f956 = 0
 flow_conservation: f176 - f716 + f376 - f736 + f476 - f746 + f576 - f756
   - f796 + f976 = 0
 flow_conservation: f196 - f916 + f396 - f936 + f496 - f946 + f596 - f956
   + f796 - f976 + f896 - f986 = 0
 flow_conservation: - f037 + f307 - f057 + f507 - f067 + f607 - f087 + f807
   = 0
 flow_conservation: - f147 + f417 - f157 + f517 - f177 + f717 - f197 + f917
   = 0
 flow_conservation: = 0
 flow_conservation: f037 - f307 - f347 + f437 - f357 + f537 - f377 + f737
   - f387 + f837 - f397 + f937 = 0
 flow_conservation: f147 - f417 + f347 - f437 - f457 + f547 - f477 + f747
   - f497 + f947 = 0
 flow_conservation: f057 - f507 + f157 - f517 + f357 - f537 + f457 - f547
   - f567 + f657 - f577 + f757 - f587 + f857 - f597 + f957 = 0
 flow_conservation: f067 - f607 + f567 - f657 - f687 + f867 = 0
 flow_conservation: f197 - f917 + f397 - f937 + f497 - f947 + f597 - f957
   + f797 - f977 + f897 - f987 = 0
 flow_conservation: - f039 + f309 - f059 + f509 - f069 + f609 - f089 + f809
   = 0
 flow_conservation: - f149 + f419 - f159 + f519 - f179 + f719 - f199 + f919
   = 0
 flow_conservation: = 0
 flow_conservation: f039 - f309 - f349 + f439 - f359 + f539 - f379 + f739
   - f389 + f839 - f399 + f939 = 0
 flow_conservation: f149 - f419 + f349 - f439 - f459 + f549 - f479 + f749
   - f499 + f949 = 0
 flow_conservation: f059 - f509 + f159 - f519 + f359 - f539 + f459 - f549
   - f569 + f659 - f579 + f759 - f589 + f859 - f599 + f959 = 0
 flow_conservation: f069 - f609 + f569 - f659 - f689 + f869 = 0
 flow_conservation: f179 - f719 + f379 - f739 + f479 - f749 + f579 - f759
   - f799 + f979 = 0
 flow_on_terminals: f141 - f411 + f151 - f511 + f171 - f711 + f191 - f911
   = -1
 flow_on_terminals: - f144 + f414 - f344 + f434 + f454 - f544 + f474 - f744
   + f494 - f944 = -1
 flow_on_terminals: - f055 + f505 - f155 + f515 - f355 + f535 - f455 + f545
   + f565 - f655 + f575 - f755 + f585 - f855 + f595 - f955 = -1
 flow_on_terminals: f030 - f300 + f050 - f500 + f060 - f600 + f080 - f800
   = -1
 flow_on_terminals: = -1
 flow_on_terminals: - f033 + f303 + f343 - f433 + f353 - f533 + f373 - f733
   + f383 - f833 + f393 - f933 = -1
 flow_on_terminals: - f066 + f606 - f566 + f656 + f686 - f866 = -1
 flow_on_terminals: - f177 + f717 - f377 + f737 - f477 + f747 - f577 + f757
   + f797 - f977 = -1
 flow_on_terminals: - f199 + f919 - f399 + f939 - f499 + f949 - f599 + f959
   - f799 + f979 - f899 + f989 = -1
 f_and_y_relation: - y03 + f031 <= 0
 f_and_y_relation: - y03 + f034 <= 0
 f_and_y_relation: - y03 + f035 <= 0
 f_and_y_relation: - y03 + f030 <= 0
 f_and_y_relation: - y03 + f032 <= 0
 f_and_y_relation: - y03 + f033 <= 0
 f_and_y_relation: - y03 + f036 <= 0
 f_and_y_relation: - y03 + f037 <= 0
 f_and_y_relation: - y03 + f039 <= 0
 f_and_y_relation: - y30 + f301 <= 0
 f_and_y_relation: - y30 + f304 <= 0
 f_and_y_relation: - y30 + f305 <= 0
 f_and_y_relation: - y30 + f300 <= 0
 f_and_y_relation: - y30 + f302 <= 0
 f_and_y_relation: - y30 + f303 <= 0
 f_and_y_relation: - y30 + f306 <= 0
 f_and_y_relation: - y30 + f307 <= 0
 f_and_y_relation: - y30 + f309 <= 0
 f_and_y_relation: - y05 + f051 <= 0
 f_and_y_relation: - y05 + f054 <= 0
 f_and_y_relation: - y05 + f055 <= 0
 f_and_y_relation: - y05 + f050 <= 0
 f_and_y_relation: - y05 + f052 <= 0
 f_and_y_relation: - y05 + f053 <= 0
 f_and_y_relation: - y05 + f056 <= 0
 f_and_y_relation: - y05 + f057 <= 0
 f_and_y_relation: - y05 + f059 <= 0
 f_and_y_relation: - y50 + f501 <= 0
 f_and_y_relation: - y50 + f504 <= 0
 f_and_y_relation: - y50 + f505 <= 0
 f_and_y_relation: - y50 + f500 <= 0
 f_and_y_relation: - y50 + f502 <= 0
 f_and_y_relation: - y50 + f503 <= 0
 f_and_y_relation: - y50 + f506 <= 0
 f_and_y_relation: - y50 + f507 <= 0
 f_and_y_relation: - y50 + f509 <= 0
 f_and_y_relation: - y06 + f061 <= 0
 f_and_y_relation: - y06 + f064 <= 0
 f_and_y_relation: - y06 + f065 <= 0
 f_and_y_relation: - y06 + f060 <= 0
 f_and_y_relation: - y06 + f062 <= 0
 f_and_y_relation: - y06 + f063 <= 0
 f_and_y_relation: - y06 + f066 <= 0
 f_and_y_relation: - y06 + f067 <= 0
 f_and_y_relation: - y06 + f069 <= 0
 f_and_y_relation: - y60 + f601 <= 0
 f_and_y_relation: - y60 + f604 <= 0
 f_and_y_relation: - y60 + f605 <= 0
 f_and_y_relation: - y60 + f600 <= 0
 f_and_y_relation: - y60 + f602 <= 0
 f_and_y_relation: - y60 + f603 <= 0
 f_and_y_relation: - y60 + f606 <= 0
 f_and_y_relation: - y60 + f607 <= 0
 f_and_y_relation: - y60 + f609 <= 0
 f_and_y_relation: - y08 + f081 <= 0
 f_and_y_relation: - y08 + f084 <= 0
 f_and_y_relation: - y08 + f085 <= 0
 f_and_y_relation: - y08 + f080 <= 0
 f_and_y_relation: - y08 + f082 <= 0
 f_and_y_relation: - y08 + f083 <= 0
 f_and_y_relation: - y08 + f086 <= 0
 f_and_y_relation: - y08 + f087 <= 0
 f_and_y_relation: - y08 + f089 <= 0
 f_and_y_relation: - y80 + f801 <= 0
 f_and_y_relation: - y80 + f804 <= 0
 f_and_y_relation: - y80 + f805 <= 0
 f_and_y_relation: - y80 + f800 <= 0
 f_and_y_relation: - y80 + f802 <= 0
 f_and_y_relation: - y80 + f803 <= 0
 f_and_y_relation: - y80 + f806 <= 0
 f_and_y_relation: - y80 + f807 <= 0
 f_and_y_relation: - y80 + f809 <= 0
 f_and_y_relation: - y14 + f141 <= 0
 f_and_y_relation: - y14 + f144 <= 0
 f_and_y_relation: - y14 + f145 <= 0
 f_and_y_relation: - y14 + f140 <= 0
 f_and_y_relation: - y14 + f142 <= 0
 f_and_y_relation: - y14 + f143 <= 0
 f_and_y_relation: - y14 + f146 <= 0
 f_and_y_relation: - y14 + f147 <= 0
 f_and_y_relation: - y14 + f149 <= 0
 f_and_y_relation: - y41 + f411 <= 0
 f_and_y_relation: - y41 + f414 <= 0
 f_and_y_relation: - y41 + f415 <= 0
 f_and_y_relation: - y41 + f410 <= 0
 f_and_y_relation: - y41 + f412 <= 0
 f_and_y_relation: - y41 + f413 <= 0
 f_and_y_relation: - y41 + f416 <= 0
 f_and_y_relation: - y41 + f417 <= 0
 f_and_y_relation: - y41 + f419 <= 0
 f_and_y_relation: - y15 + f151 <= 0
 f_and_y_relation: - y15 + f154 <= 0
 f_and_y_relation: - y15 + f155 <= 0
 f_and_y_relation: - y15 + f150 <= 0
 f_and_y_relation: - y15 + f152 <= 0
 f_and_y_relation: - y15 + f153 <= 0
 f_and_y_relation: - y15 + f156 <= 0
 f_and_y_relation: - y15 + f157 <= 0
 f_and_y_relation: - y15 + f159 <= 0
 f_and_y_relation: - y51 + f511 <= 0
 f_and_y_relation: - y51 + f514 <= 0
 f_and_y_relation: - y51 + f515 <= 0
 f_and_y_relation: - y51 + f510 <= 0
 f_and_y_relation: - y51 + f512 <= 0
 f_and_y_relation: - y51 + f513 <= 0
 f_and_y_relation: - y51 + f516 <= 0
 f_and_y_relation: - y51 + f517 <= 0
 f_and_y_relation: - y51 + f519 <= 0
 f_and_y_relation: - y17 + f171 <= 0
 f_and_y_relation: - y17 + f174 <= 0
 f_and_y_relation: - y17 + f175 <= 0
 f_and_y_relation: - y17 + f170 <= 0
 f_and_y_relation: - y17 + f172 <= 0
 f_and_y_relation: - y17 + f173 <= 0
 f_and_y_relation: - y17 + f176 <= 0
 f_and_y_relation: - y17 + f177 <= 0
 f_and_y_relation: - y17 + f179 <= 0
 f_and_y_relation: - y71 + f711 <= 0
 f_and_y_relation: - y71 + f714 <= 0
 f_and_y_relation: - y71 + f715 <= 0
 f_and_y_relation: - y71 + f710 <= 0
 f_and_y_relation: - y71 + f712 <= 0
 f_and_y_relation: - y71 + f713 <= 0
 f_and_y_relation: - y71 + f716 <= 0
 f_and_y_relation: - y71 + f717 <= 0
 f_and_y_relation: - y71 + f719 <= 0
 f_and_y_relation: - y19 + f191 <= 0
 f_and_y_relation: - y19 + f194 <= 0
 f_and_y_relation: - y19 + f195 <= 0
 f_and_y_relation: - y19 + f190 <= 0
 f_and_y_relation: - y19 + f192 <= 0
 f_and_y_relation: - y19 + f193 <= 0
 f_and_y_relation: - y19 + f196 <= 0
 f_and_y_relation: - y19 + f197 <= 0
 f_and_y_relation: - y19 + f199 <= 0
 f_and_y_relation: - y91 + f911 <= 0
 f_and_y_relation: - y91 + f914 <= 0
 f_and_y_relation: - y91 + f915 <= 0
 f_and_y_relation: - y91 + f910 <= 0
 f_and_y_relation: - y91 + f912 <= 0
 f_and_y_relation: - y91 + f913 <= 0
 f_and_y_relation: - y91 + f916 <= 0
 f_and_y_relation: - y91 + f917 <= 0
 f_and_y_relation: - y91 + f919 <= 0
 f_and_y_relation: - y34 + f341 <= 0
 f_and_y_relation: - y34 + f344 <= 0
 f_and_y_relation: - y34 + f345 <= 0
 f_and_y_relation: - y34 + f340 <= 0
 f_and_y_relation: - y34 + f342 <= 0
 f_and_y_relation: - y34 + f343 <= 0
 f_and_y_relation: - y34 + f346 <= 0
 f_and_y_relation: - y34 + f347 <= 0
 f_and_y_relation: - y34 + f349 <= 0
 f_and_y_relation: - y43 + f431 <= 0
 f_and_y_relation: - y43 + f434 <= 0
 f_and_y_relation: - y43 + f435 <= 0
 f_and_y_relation: - y43 + f430 <= 0
 f_and_y_relation: - y43 + f432 <= 0
 f_and_y_relation: - y43 + f433 <= 0
 f_and_y_relation: - y43 + f436 <= 0
 f_and_y_relation: - y43 + f437 <= 0
 f_and_y_relation: - y43 + f439 <= 0
 f_and_y_relation: - y35 + f351 <= 0
 f_and_y_relation: - y35 + f354 <= 0
 f_and_y_relation: - y35 + f355 <= 0
 f_and_y_relation: - y35 + f350 <= 0
 f_and_y_relation: - y35 + f352 <= 0
 f_and_y_relation: - y35 + f353 <= 0
 f_and_y_relation: - y35 + f356 <= 0
 f_and_y_relation: - y35 + f357 <= 0
 f_and_y_relation: - y35 + f359 <= 0
 f_and_y_relation: - y53 + f531 <= 0
 f_and_y_relation: - y53 + f534 <= 0
 f_and_y_relation: - y53 + f535 <= 0
 f_and_y_relation: - y53 + f530 <= 0
 f_and_y_relation: - y53 + f532 <= 0
 f_and_y_relation: - y53 + f533 <= 0
 f_and_y_relation: - y53 + f536 <= 0
 f_and_y_relation: - y53 + f537 <= 0
 f_and_y_relation: - y53 + f539 <= 0
 f_and_y_relation: - y37 + f371 <= 0
 f_and_y_relation: - y37 + f374 <= 0
 f_and_y_relation: - y37 + f375 <= 0
 f_and_y_relation: - y37 + f370 <= 0
 f_and_y_relation: - y37 + f372 <= 0
 f_and_y_relation: - y37 + f373 <= 0
 f_and_y_relation: - y37 + f376 <= 0
 f_and_y_relation: - y37 + f377 <= 0
 f_and_y_relation: - y37 + f379 <= 0
 f_and_y_relation: - y73 + f731 <= 0
 f_and_y_relation: - y73 + f734 <= 0
 f_and_y_relation: - y73 + f735 <= 0
 f_and_y_relation: - y73 + f730 <= 0
 f_and_y_relation: - y73 + f732 <= 0
 f_and_y_relation: - y73 + f733 <= 0
 f_and_y_relation: - y73 + f736 <= 0
 f_and_y_relation: - y73 + f737 <= 0
 f_and_y_relation: - y73 + f739 <= 0
 f_and_y_relation: - y38 + f381 <= 0
 f_and_y_relation: - y38 + f384 <= 0
 f_and_y_relation: - y38 + f385 <= 0
 f_and_y_relation: - y38 + f380 <= 0
 f_and_y_relation: - y38 + f382 <= 0
 f_and_y_relation: - y38 + f383 <= 0
 f_and_y_relation: - y38 + f386 <= 0
 f_and_y_relation: - y38 + f387 <= 0
 f_and_y_relation: - y38 + f389 <= 0
 f_and_y_relation: - y83 + f831 <= 0
 f_and_y_relation: - y83 + f834 <= 0
 f_and_y_relation: - y83 + f835 <= 0
 f_and_y_relation: - y83 + f830 <= 0
 f_and_y_relation: - y83 + f832 <= 0
 f_and_y_relation: - y83 + f833 <= 0
 f_and_y_relation: - y83 + f836 <= 0
 f_and_y_relation: - y83 + f837 <= 0
 f_and_y_relation: - y83 + f839 <= 0
 f_and_y_relation: - y39 + f391 <= 0
 f_and_y_relation: - y39 + f394 <= 0
 f_and_y_relation: - y39 + f395 <= 0
 f_and_y_relation: - y39 + f390 <= 0
 f_and_y_relation: - y39 + f392 <= 0
 f_and_y_relation: - y39 + f393 <= 0
 f_and_y_relation: - y39 + f396 <= 0
 f_and_y_relation: - y39 + f397 <= 0
 f_and_y_relation: - y39 + f399 <= 0
 f_and_y_relation: - y93 + f931 <= 0
 f_and_y_relation: - y93 + f934 <= 0
 f_and_y_relation: - y93 + f935 <= 0
 f_and_y_relation: - y93 + f930 <= 0
 f_and_y_relation: - y93 + f932 <= 0
 f_and_y_relation: - y93 + f933 <= 0
 f_and_y_relation: - y93 + f936 <= 0
 f_and_y_relation: - y93 + f937 <= 0
 f_and_y_relation: - y93 + f939 <= 0
 f_and_y_relation: - y45 + f451 <= 0
 f_and_y_relation: - y45 + f454 <= 0
 f_and_y_relation: - y45 + f455 <= 0
 f_and_y_relation: - y45 + f450 <= 0
 f_and_y_relation: - y45 + f452 <= 0
 f_and_y_relation: - y45 + f453 <= 0
 f_and_y_relation: - y45 + f456 <= 0
 f_and_y_relation: - y45 + f457 <= 0
 f_and_y_relation: - y45 + f459 <= 0
 f_and_y_relation: - y54 + f541 <= 0
 f_and_y_relation: - y54 + f544 <= 0
 f_and_y_relation: - y54 + f545 <= 0
 f_and_y_relation: - y54 + f540 <= 0
 f_and_y_relation: - y54 + f542 <= 0
 f_and_y_relation: - y54 + f543 <= 0
 f_and_y_relation: - y54 + f546 <= 0
 f_and_y_relation: - y54 + f547 <= 0
 f_and_y_relation: - y54 + f549 <= 0
 f_and_y_relation: - y47 + f471 <= 0
 f_and_y_relation: - y47 + f474 <= 0
 f_and_y_relation: - y47 + f475 <= 0
 f_and_y_relation: - y47 + f470 <= 0
 f_and_y_relation: - y47 + f472 <= 0
 f_and_y_relation: - y47 + f473 <= 0
 f_and_y_relation: - y47 + f476 <= 0
 f_and_y_relation: - y47 + f477 <= 0
 f_and_y_relation: - y47 + f479 <= 0
 f_and_y_relation: - y74 + f741 <= 0
 f_and_y_relation: - y74 + f744 <= 0
 f_and_y_relation: - y74 + f745 <= 0
 f_and_y_relation: - y74 + f740 <= 0
 f_and_y_relation: - y74 + f742 <= 0
 f_and_y_relation: - y74 + f743 <= 0
 f_and_y_relation: - y74 + f746 <= 0
 f_and_y_relation: - y74 + f747 <= 0
 f_and_y_relation: - y74 + f749 <= 0
 f_and_y_relation: - y49 + f491 <= 0
 f_and_y_relation: - y49 + f494 <= 0
 f_and_y_relation: - y49 + f495 <= 0
 f_and_y_relation: - y49 + f490 <= 0
 f_and_y_relation: - y49 + f492 <= 0
 f_and_y_relation: - y49 + f493 <= 0
 f_and_y_relation: - y49 + f496 <= 0
 f_and_y_relation: - y49 + f497 <= 0
 f_and_y_relation: - y49 + f499 <= 0
 f_and_y_relation: - y94 + f941 <= 0
 f_and_y_relation: - y94 + f944 <= 0
 f_and_y_relation: - y94 + f945 <= 0
 f_and_y_relation: - y94 + f940 <= 0
 f_and_y_relation: - y94 + f942 <= 0
 f_and_y_relation: - y94 + f943 <= 0
 f_and_y_relation: - y94 + f946 <= 0
 f_and_y_relation: - y94 + f947 <= 0
 f_and_y_relation: - y94 + f949 <= 0
 f_and_y_relation: - y56 + f561 <= 0
 f_and_y_relation: - y56 + f564 <= 0
 f_and_y_relation: - y56 + f565 <= 0
 f_and_y_relation: - y56 + f560 <= 0
 f_and_y_relation: - y56 + f562 <= 0
 f_and_y_relation: - y56 + f563 <= 0
 f_and_y_relation: - y56 + f566 <= 0
 f_and_y_relation: - y56 + f567 <= 0
 f_and_y_relation: - y56 + f569 <= 0
 f_and_y_relation: - y65 + f651 <= 0
 f_and_y_relation: - y65 + f654 <= 0
 f_and_y_relation: - y65 + f655 <= 0
 f_and_y_relation: - y65 + f650 <= 0
 f_and_y_relation: - y65 + f652 <= 0
 f_and_y_relation: - y65 + f653 <= 0
 f_and_y_relation: - y65 + f656 <= 0
 f_and_y_relation: - y65 + f657 <= 0
 f_and_y_relation: - y65 + f659 <= 0
 f_and_y_relation: - y57 + f571 <= 0
 f_and_y_relation: - y57 + f574 <= 0
 f_and_y_relation: - y57 + f575 <= 0
 f_and_y_relation: - y57 + f570 <= 0
 f_and_y_relation: - y57 + f572 <= 0
 f_and_y_relation: - y57 + f573 <= 0
 f_and_y_relation: - y57 + f576 <= 0
 f_and_y_relation: - y57 + f577 <= 0
 f_and_y_relation: - y57 + f579 <= 0
 f_and_y_relation: - y75 + f751 <= 0
 f_and_y_relation: - y75 + f754 <= 0
 f_and_y_relation: - y75 + f755 <= 0
 f_and_y_relation: - y75 + f750 <= 0
 f_and_y_relation: - y75 + f752 <= 0
 f_and_y_relation: - y75 + f753 <= 0
 f_and_y_relation: - y75 + f756 <= 0
 f_and_y_relation: - y75 + f757 <= 0
 f_and_y_relation: - y75 + f759 <= 0
 f_and_y_relation: - y58 + f581 <= 0
 f_and_y_relation: - y58 + f584 <= 0
 f_and_y_relation: - y58 + f585 <= 0
 f_and_y_relation: - y58 + f580 <= 0
 f_and_y_relation: - y58 + f582 <= 0
 f_and_y_relation: - y58 + f583 <= 0
 f_and_y_relation: - y58 + f586 <= 0
 f_and_y_relation: - y58 + f587 <= 0
 f_and_y_relation: - y58 + f589 <= 0
 f_and_y_relation: - y85 + f851 <= 0
 f_and_y_relation: - y85 + f854 <= 0
 f_and_y_relation: - y85 + f855 <= 0
 f_and_y_relation: - y85 + f850 <= 0
 f_and_y_relation: - y85 + f852 <= 0
 f_and_y_relation: - y85 + f853 <= 0
 f_and_y_relation: - y85 + f856 <= 0
 f_and_y_relation: - y85 + f857 <= 0
 f_and_y_relation: - y85 + f859 <= 0
 f_and_y_relation: - y59 + f591 <= 0
 f_and_y_relation: - y59 + f594 <= 0
 f_and_y_relation: - y59 + f595 <= 0
 f_and_y_relation: - y59 + f590 <= 0
 f_and_y_relation: - y59 + f592 <= 0
 f_and_y_relation: - y59 + f593 <= 0
 f_and_y_relation: - y59 + f596 <= 0
 f_and_y_relation: - y59 + f597 <= 0
 f_and_y_relation: - y59 + f599 <= 0
 f_and_y_relation: - y95 + f951 <= 0
 f_and_y_relation: - y95 + f954 <= 0
 f_and_y_relation: - y95 + f955 <= 0
 f_and_y_relation: - y95 + f950 <= 0
 f_and_y_relation: - y95 + f952 <= 0
 f_and_y_relation: - y95 + f953 <= 0
 f_and_y_relation: - y95 + f956 <= 0
 f_and_y_relation: - y95 + f957 <= 0
 f_and_y_relation: - y95 + f959 <= 0
 f_and_y_relation: - y68 + f681 <= 0
 f_and_y_relation: - y68 + f684 <= 0
 f_and_y_relation: - y68 + f685 <= 0
 f_and_y_relation: - y68 + f680 <= 0
 f_and_y_relation: - y68 + f682 <= 0
 f_and_y_relation: - y68 + f683 <= 0
 f_and_y_relation: - y68 + f686 <= 0
 f_and_y_relation: - y68 + f687 <= 0
 f_and_y_relation: - y68 + f689 <= 0
 f_and_y_relation: - y86 + f861 <= 0
 f_and_y_relation: - y86 + f864 <= 0
 f_and_y_relation: - y86 + f865 <= 0
 f_and_y_relation: - y86 + f860 <= 0
 f_and_y_relation: - y86 + f862 <= 0
 f_and_y_relation: - y86 + f863 <= 0
 f_and_y_relation: - y86 + f866 <= 0
 f_and_y_relation: - y86 + f867 <= 0
 f_and_y_relation: - y86 + f869 <= 0
 f_and_y_relation: - y79 + f791 <= 0
 f_and_y_relation: - y79 + f794 <= 0
 f_and_y_relation: - y79 + f795 <= 0
 f_and_y_relation: - y79 + f790 <= 0
 f_and_y_relation: - y79 + f792 <= 0
 f_and_y_relation: - y79 + f793 <= 0
 f_and_y_relation: - y79 + f796 <= 0
 f_and_y_relation: - y79 + f797 <= 0
 f_and_y_relation: - y79 + f799 <= 0
 f_and_y_relation: - y97 + f971 <= 0
 f_and_y_relation: - y97 + f974 <= 0
 f_and_y_relation: - y97 + f975 <= 0
 f_and_y_relation: - y97 + f970 <= 0
 f_and_y_relation: - y97 + f972 <= 0
 f_and_y_relation: - y97 + f973 <= 0
 f_and_y_relation: - y97 + f976 <= 0
 f_and_y_relation: - y97 + f977 <= 0
 f_and_y_relation: - y97 + f979 <= 0
 f_and_y_relation: - y89 + f891 <= 0
 f_and_y_relation: - y89 + f894 <= 0
 f_and_y_relation: - y89 + f895 <= 0
 f_and_y_relation: - y89 + f890 <= 0
 f_and_y_relation: - y89 + f892 <= 0
 f_and_y_relation: - y89 + f893 <= 0
 f_and_y_relation: - y89 + f896 <= 0
 f_and_y_relation: - y89 + f897 <= 0
 f_and_y_relation: - y89 + f899 <= 0
 f_and_y_relation: - y98 + f981 <= 0
 f_and_y_relation: - y98 + f984 <= 0
 f_and_y_relation: - y98 + f985 <= 0
 f_and_y_relation: - y98 + f980 <= 0
 f_and_y_relation: - y98 + f982 <= 0
 f_and_y_relation: - y98 + f983 <= 0
 f_and_y_relation: - y98 + f986 <= 0
 f_and_y_relation: - y98 + f987 <= 0
 f_and_y_relation: - y98 + f989 <= 0
 maximum_of_arcs: y03 + y30 + y05 + y50 + y06 + y60 + y08 + y80 + y14 + y41
   + y15 + y51 + y17 + y71 + y19 + y91 + y34 + y43 + y35 + y53 + y37 + y73
   + y38 + y83 + y39 + y93 + y45 + y54 + y47 + y74 + y49 + y94 + y56 + y65
   + y57 + y75 + y58 + y85 + y59 + y95 + y68 + y86 + y79 + y97 + y89 + y98
   = 9
 delay_limit: 2767 f031 + 2767 f301 + 3414 f051 + 3414 f501 + 2594 f061
   + 2594 f601 + 2866 f081 + 2866 f801 + 2593 f141 + 2593 f411 + 3451 f151
   + 3451 f511 + 2565 f171 + 2565 f711 + 2698 f191 + 2698 f911 + 2946 f341
   + 2946 f431 + 3404 f351 + 3404 f531 + 2949 f371 + 2949 f731 + 2974 f381
   + 2974 f831 + 2820 f391 + 2820 f931 + 3379 f451 + 3379 f541 + 2454 f471
   + 2454 f741 + 2721 f491 + 2721 f941 + 3658 f561 + 3658 f651 + 3387 f571
   + 3387 f751 + 3404 f581 + 3404 f851 + 3404 f591 + 3404 f951 + 3093 f681
   + 3093 f861 + 2727 f791 + 2727 f971 + 2862 f891 + 2862 f981 - 18173 z1
   <= 5766
 jitter_limit: 3974 f031 + 3974 f301 + 5682 f051 + 5682 f501 + 4983 f061
   + 4983 f601 + 3870 f081 + 3870 f801 + 2689 f141 + 2689 f411 + 5578 f151
   + 5578 f511 + 2714 f171 + 2714 f711 + 2940 f191 + 2940 f911 + 763 f341
   + 763 f431 + 5105 f351 + 5105 f531 + 340 f371 + 340 f731 + 3100 f381
   + 3100 f831 + 1326 f391 + 1326 f931 + 5374 f451 + 5374 f541 + 927 f471
   + 927 f741 + 897 f491 + 897 f941 + 4560 f561 + 4560 f651 + 5672 f571
   + 5672 f751 + 5105 f581 + 5105 f851 + 5105 f591 + 5105 f951 + 1357 f681
   + 1357 f861 + 1044 f791 + 1044 f971 + 1490 f891 + 1490 f981 - 6019 z1
   <= 4814
 delay_limit: 2767 f034 + 2767 f304 + 3414 f054 + 3414 f504 + 2594 f064
   + 2594 f604 + 2866 f084 + 2866 f804 + 2593 f144 + 2593 f414 + 3451 f154
   + 3451 f514 + 2565 f174 + 2565 f714 + 2698 f194 + 2698 f914 + 2946 f344
   + 2946 f434 + 3404 f354 + 3404 f534 + 2949 f374 + 2949 f734 + 2974 f384
   + 2974 f834 + 2820 f394 + 2820 f934 + 3379 f454 + 3379 f544 + 2454 f474
   + 2454 f744 + 2721 f494 + 2721 f944 + 3658 f564 + 3658 f654 + 3387 f574
   + 3387 f754 + 3404 f584 + 3404 f854 + 3404 f594 + 3404 f954 + 3093 f684
   + 3093 f864 + 2727 f794 + 2727 f974 + 2862 f894 + 2862 f984 - 18173 z4
   <= 5766
 jitter_limit: 3974 f034 + 3974 f304 + 5682 f054 + 5682 f504 + 4983 f064
   + 4983 f604 + 3870 f084 + 3870 f804 + 2689 f144 + 2689 f414 + 5578 f154
   + 5578 f514 + 2714 f174 + 2714 f714 + 2940 f194 + 2940 f914 + 763 f344
   + 763 f434 + 5105 f354 + 5105 f534 + 340 f374 + 340 f734 + 3100 f384
   + 3100 f834 + 1326 f394 + 1326 f934 + 5374 f454 + 5374 f544 + 927 f474
   + 927 f744 + 897 f494 + 897 f944 + 4560 f564 + 4560 f654 + 5672 f574
   + 5672 f754 + 5105 f584 + 5105 f854 + 5105 f594 + 5105 f954 + 1357 f684
   + 1357 f864 + 1044 f794 + 1044 f974 + 1490 f894 + 1490 f984 - 6019 z4
   <= 4814
 delay_limit: 2767 f035 + 2767 f305 + 3414 f055 + 3414 f505 + 2594 f065
   + 2594 f605 + 2866 f085 + 2866 f805 + 2593 f145 + 2593 f415 + 3451 f155
   + 3451 f515 + 2565 f175 + 2565 f715 + 2698 f195 + 2698 f915 + 2946 f345
   + 2946 f435 + 3404 f355 + 3404 f535 + 2949 f375 + 2949 f735 + 2974 f385
   + 2974 f835 + 2820 f395 + 2820 f935 + 3379 f455 + 3379 f545 + 2454 f475
   + 2454 f745 + 2721 f495 + 2721 f945 + 3658 f565 + 3658 f655 + 3387 f575
   + 3387 f755 + 3404 f585 + 3404 f855 + 3404 f595 + 3404 f955 + 3093 f685
   + 3093 f865 + 2727 f795 + 2727 f975 + 2862 f895 + 2862 f985 - 18173 z5
   <= 5766
 jitter_limit: 3974 f035 + 3974 f305 + 5682 f055 + 5682 f505 + 4983 f065
   + 4983 f605 + 3870 f085 + 3870 f805 + 2689 f145 + 2689 f415 + 5578 f155
   + 5578 f515 + 2714 f175 + 2714 f715 + 2940 f195 + 2940 f915 + 763 f345
   + 763 f435 + 5105 f355 + 5105 f535 + 340 f375 + 340 f735 + 3100 f385
   + 3100 f835 + 1326 f395 + 1326 f935 + 5374 f455 + 5374 f545 + 927 f475
   + 927 f745 + 897 f495 + 897 f945 + 4560 f565 + 4560 f655 + 5672 f575
   + 5672 f755 + 5105 f585 + 5105 f855 + 5105 f595 + 5105 f955 + 1357 f685
   + 1357 f865 + 1044 f795 + 1044 f975 + 1490 f895 + 1490 f985 - 6019 z5
   <= 4814
 limit_of_variation_between_pairs: 2767 f031 - 2767 f034 + 2767 f301
   - 2767 f304 + 3414 f051 - 3414 f054 + 3414 f501 - 3414 f504 + 2594 f061
   - 2594 f064 + 2594 f601 - 2594 f604 + 2866 f081 - 2866 f084 + 2866 f801
   - 2866 f804 + 2593 f141 - 2593 f144 + 2593 f411 - 2593 f414 + 3451 f151
   - 3451 f154 + 3451 f511 - 3451 f514 + 2565 f171 - 2565 f174 + 2565 f711
   - 2565 f714 + 2698 f191 - 2698 f194 + 2698 f911 - 2698 f914 + 2946 f341
   - 2946 f344 + 2946 f431 - 2946 f434 + 3404 f351 - 3404 f354 + 3404 f531
   - 3404 f534 + 2949 f371 - 2949 f374 + 2949 f731 - 2949 f734 + 2974 f381
   - 2974 f384 + 2974 f831 - 2974 f834 + 2820 f391 - 2820 f394 + 2820 f931
   - 2820 f934 + 3379 f451 - 3379 f454 + 3379 f541 - 3379 f544 + 2454 f471
   - 2454 f474 + 2454 f741 - 2454 f744 + 2721 f491 - 2721 f494 + 2721 f941
   - 2721 f944 + 3658 f561 - 3658 f564 + 3658 f651 - 3658 f654 + 3387 f571
   - 3387 f574 + 3387 f751 - 3387 f754 + 3404 f581 - 3404 f584 + 3404 f851
   - 3404 f854 + 3404 f591 - 3404 f594 + 3404 f951 - 3404 f954 + 3093 f681
   - 3093 f684 + 3093 f861 - 3093 f864 + 2727 f791 - 2727 f794 + 2727 f971
   - 2727 f974 + 2862 f891 - 2862 f894 + 2862 f981 - 2862 f984 - 18173 z1
   + 3128 z4 <= 3311
 limit_of_variation_between_pairs: 2767 f031 - 2767 f035 + 2767 f301
   - 2767 f305 + 3414 f051 - 3414 f055 + 3414 f501 - 3414 f505 + 2594 f061
   - 2594 f065 + 2594 f601 - 2594 f605 + 2866 f081 - 2866 f085 + 2866 f801
   - 2866 f805 + 2593 f141 - 2593 f145 + 2593 f411 - 2593 f415 + 3451 f151
   - 3451 f155 + 3451 f511 - 3451 f515 + 2565 f171 - 2565 f175 + 2565 f711
   - 2565 f715 + 2698 f191 - 2698 f195 + 2698 f911 - 2698 f915 + 2946 f341
   - 2946 f345 + 2946 f431 - 2946 f435 + 3404 f351 - 3404 f355 + 3404 f531
   - 3404 f535 + 2949 f371 - 2949 f375 + 2949 f731 - 2949 f735 + 2974 f381
   - 2974 f385 + 2974 f831 - 2974 f835 + 2820 f391 - 2820 f395 + 2820 f931
   - 2820 f935 + 3379 f451 - 3379 f455 + 3379 f541 - 3379 f545 + 2454 f471
   - 2454 f475 + 2454 f741 - 2454 f745 + 2721 f491 - 2721 f495 + 2721 f941
   - 2721 f945 + 3658 f561 - 3658 f565 + 3658 f651 - 3658 f655 + 3387 f571
   - 3387 f575 + 3387 f751 - 3387 f755 + 3404 f581 - 3404 f585 + 3404 f851
   - 3404 f855 + 3404 f591 - 3404 f595 + 3404 f951 - 3404 f955 + 3093 f681
   - 3093 f685 + 3093 f861 - 3093 f865 + 2727 f791 - 2727 f795 + 2727 f971
   - 2727 f975 + 2862 f891 - 2862 f895 + 2862 f981 - 2862 f985 - 18173 z1
   + 949 z5 <= 3311
 limit_of_variation_between_pairs: - 2767 f031 + 2767 f034 - 2767 f301
   + 2767 f304 - 3414 f051 + 3414 f054 - 3414 f501 + 3414 f504 - 2594 f061
   + 2594 f064 - 2594 f601 + 2594 f604 - 2866 f081 + 2866 f084 - 2866 f801
   + 2866 f804 - 2593 f141 + 2593 f144 - 2593 f411 + 2593 f414 - 3451 f151
   + 3451 f154 - 3451 f511 + 3451 f514 - 2565 f171 + 2565 f174 - 2565 f711
   + 2565 f714 - 2698 f191 + 2698 f194 - 2698 f911 + 2698 f914 - 2946 f341
   + 2946 f344 - 2946 f431 + 2946 f434 - 3404 f351 + 3404 f354 - 3404 f531
   + 3404 f534 - 2949 f371 + 2949 f374 - 2949 f731 + 2949 f734 - 2974 f381
   + 2974 f384 - 2974 f831 + 2974 f834 - 2820 f391 + 2820 f394 - 2820 f931
   + 2820 f934 - 3379 f451 + 3379 f454 - 3379 f541 + 3379 f544 - 2454 f471
   + 2454 f474 - 2454 f741 + 2454 f744 - 2721 f491 + 2721 f494 - 2721 f941
   + 2721 f944 - 3658 f561 + 3658 f564 - 3658 f651 + 3658 f654 - 3387 f571
   + 3387 f574 - 3387 f751 + 3387 f754 - 3404 f581 + 3404 f584 - 3404 f851
   + 3404 f854 - 3404 f591 + 3404 f594 - 3404 f951 + 3404 f954 - 3093 f681
   + 3093 f684 - 3093 f861 + 3093 f864 - 2727 f791 + 2727 f794 - 2727 f971
   + 2727 f974 - 2862 f891 + 2862 f894 - 2862 f981 + 2862 f984 + 3105 z1
   - 18173 z4 <= 3311
 limit_of_variation_between_pairs: 2767 f034 - 2767 f035 + 2767 f304
   - 2767 f305 + 3414 f054 - 3414 f055 + 3414 f504 - 3414 f505 + 2594 f064
   - 2594 f065 + 2594 f604 - 2594 f605 + 2866 f084 - 2866 f085 + 2866 f804
   - 2866 f805 + 2593 f144 - 2593 f145 + 2593 f414 - 2593 f415 + 3451 f154
   - 3451 f155 + 3451 f514 - 3451 f515 + 2565 f174 - 2565 f175 + 2565 f714
   - 2565 f715 + 2698 f194 - 2698 f195 + 2698 f914 - 2698 f915 + 2946 f344
   - 2946 f345 + 2946 f434 - 2946 f435 + 3404 f354 - 3404 f355 + 3404 f534
   - 3404 f535 + 2949 f374 - 2949 f375 + 2949 f734 - 2949 f735 + 2974 f384
   - 2974 f385 + 2974 f834 - 2974 f835 + 2820 f394 - 2820 f395 + 2820 f934
   - 2820 f935 + 3379 f454 - 3379 f455 + 3379 f544 - 3379 f545 + 2454 f474
   - 2454 f475 + 2454 f744 - 2454 f745 + 2721 f494 - 2721 f495 + 2721 f944
   - 2721 f945 + 3658 f564 - 3658 f565 + 3658 f654 - 3658 f655 + 3387 f574
   - 3387 f575 + 3387 f754 - 3387 f755 + 3404 f584 - 3404 f585 + 3404 f854
   - 3404 f855 + 3404 f594 - 3404 f595 + 3404 f954 - 3404 f955 + 3093 f684
   - 3093 f685 + 3093 f864 - 3093 f865 + 2727 f794 - 2727 f795 + 2727 f974
   - 2727 f975 + 2862 f894 - 2862 f895 + 2862 f984 - 2862 f985 - 18173 z4
   + 949 z5 <= 3311
 limit_of_variation_between_pairs: - 2767 f031 + 2767 f035 - 2767 f301
   + 2767 f305 - 3414 f051 + 3414 f055 - 3414 f501 + 3414 f505 - 2594 f061
   + 2594 f065 - 2594 f601 + 2594 f605 - 2866 f081 + 2866 f085 - 2866 f801
   + 2866 f805 - 2593 f141 + 2593 f145 - 2593 f411 + 2593 f415 - 3451 f151
   + 3451 f155 - 3451 f511 + 3451 f515 - 2565 f171 + 2565 f175 - 2565 f711
   + 2565 f715 - 2698 f191 + 2698 f195 - 2698 f911 + 2698 f915 - 2946 f341
   + 2946 f345 - 2946 f431 + 2946 f435 - 3404 f351 + 3404 f355 - 3404 f531
   + 3404 f535 - 2949 f371 + 2949 f375 - 2949 f731 + 2949 f735 - 2974 f381
   + 2974 f385 - 2974 f831 + 2974 f835 - 2820 f391 + 2820 f395 - 2820 f931
   + 2820 f935 - 3379 f451 + 3379 f455 - 3379 f541 + 3379 f545 - 2454 f471
   + 2454 f475 - 2454 f741 + 2454 f745 - 2721 f491 + 2721 f495 - 2721 f941
   + 2721 f945 - 3658 f561 + 3658 f565 - 3658 f651 + 3658 f655 - 3387 f571
   + 3387 f575 - 3387 f751 + 3387 f755 - 3404 f581 + 3404 f585 - 3404 f851
   + 3404 f855 - 3404 f591 + 3404 f595 - 3404 f951 + 3404 f955 - 3093 f681
   + 3093 f685 - 3093 f861 + 3093 f865 - 2727 f791 + 2727 f795 - 2727 f971
   + 2727 f975 - 2862 f891 + 2862 f895 - 2862 f981 + 2862 f985 + 3105 z1
   - 18173 z5 <= 3311
 limit_of_variation_between_pairs: - 2767 f034 + 2767 f035 - 2767 f304
   + 2767 f305 - 3414 f054 + 3414 f055 - 3414 f504 + 3414 f505 - 2594 f064
   + 2594 f065 - 2594 f604 + 2594 f605 - 2866 f084 + 2866 f085 - 2866 f804
   + 2866 f805 - 2593 f144 + 2593 f145 - 2593 f414 + 2593 f415 - 3451 f154
   + 3451 f155 - 3451 f514 + 3451 f515 - 2565 f174 + 2565 f175 - 2565 f714
   + 2565 f715 - 2698 f194 + 2698 f195 - 2698 f914 + 2698 f915 - 2946 f344
   + 2946 f345 - 2946 f434 + 2946 f435 - 3404 f354 + 3404 f355 - 3404 f534
   + 3404 f535 - 2949 f374 + 2949 f375 - 2949 f734 + 2949 f735 - 2974 f384
   + 2974 f385 - 2974 f834 + 2974 f835 - 2820 f394 + 2820 f395 - 2820 f934
   + 2820 f935 - 3379 f454 + 3379 f455 - 3379 f544 + 3379 f545 - 2454 f474
   + 2454 f475 - 2454 f744 + 2454 f745 - 2721 f494 + 2721 f495 - 2721 f944
   + 2721 f945 - 3658 f564 + 3658 f565 - 3658 f654 + 3658 f655 - 3387 f574
   + 3387 f575 - 3387 f754 + 3387 f755 - 3404 f584 + 3404 f585 - 3404 f854
   + 3404 f855 - 3404 f594 + 3404 f595 - 3404 f954 + 3404 f955 - 3093 f684
   + 3093 f685 - 3093 f864 + 3093 f865 - 2727 f794 + 2727 f795 - 2727 f974
   + 2727 f975 - 2862 f894 + 2862 f895 - 2862 f984 + 2862 f985 + 3128 z4
   - 18173 z5 <= 3311
Bounds
 f031 <= 1
 f034 <= 1
 f035 <= 1
 f030 <= 1
 f032 <= 1
 f033 <= 1
 f036 <= 1
 f037 <= 1
 f039 <= 1
 f301 <= 1
 f304 <= 1
 f305 <= 1
 f300 <= 1
 f302 <= 1
 f303 <= 1
 f306 <= 1
 f307 <= 1
 f309 <= 1
 f051 <= 1
 f054 <= 1
 f055 <= 1
 f050 <= 1
 f052 <= 1
 f053 <= 1
 f056 <= 1
 f057 <= 1
 f059 <= 1
 f501 <= 1
 f504 <= 1
 f505 <= 1
 f500 <= 1
 f502 <= 1
 f503 <= 1
 f506 <= 1
 f507 <= 1
 f509 <= 1
 f061 <= 1
 f064 <= 1
 f065 <= 1
 f060 <= 1
 f062 <= 1
 f063 <= 1
 f066 <= 1
 f067 <= 1
 f069 <= 1
 f601 <= 1
 f604 <= 1
 f605 <= 1
 f600 <= 1
 f602 <= 1
 f603 <= 1
 f606 <= 1
 f607 <= 1
 f609 <= 1
 f081 <= 1
 f084 <= 1
 f085 <= 1
 f080 <= 1
 f082 <= 1
 f083 <= 1
 f086 <= 1
 f087 <= 1
 f089 <= 1
 f801 <= 1
 f804 <= 1
 f805 <= 1
 f800 <= 1
 f802 <= 1
 f803 <= 1
 f806 <= 1
 f807 <= 1
 f809 <= 1
 f141 <= 1
 f144 <= 1
 f145 <= 1
 f140 <= 1
 f142 <= 1
 f143 <= 1
 f146 <= 1
 f147 <= 1
 f149 <= 1
 f411 <= 1
 f414 <= 1
 f415 <= 1
 f410 <= 1
 f412 <= 1
 f413 <= 1
 f416 <= 1
 f417 <= 1
 f419 <= 1
 f151 <= 1
 f154 <= 1
 f155 <= 1
 f150 <= 1
 f152 <= 1
 f153 <= 1
 f156 <= 1
 f157 <= 1
 f159 <= 1
 f511 <= 1
 f514 <= 1
 f515 <= 1
 f510 <= 1
 f512 <= 1
 f513 <= 1
 f516 <= 1
 f517 <= 1
 f519 <= 1
 f171 <= 1
 f174 <= 1
 f175 <= 1
 f170 <= 1
 f172 <= 1
 f173 <= 1
 f176 <= 1
 f177 <= 1
 f179 <= 1
 f711 <= 1
 f714 <= 1
 f715 <= 1
 f710 <= 1
 f712 <= 1
 f713 <= 1
 f716 <= 1
 f717 <= 1
 f719 <= 1
 f191 <= 1
 f194 <= 1
 f195 <= 1
 f190 <= 1
 f192 <= 1
 f193 <= 1
 f196 <= 1
 f197 <= 1
 f199 <= 1
 f911 <= 1
 f914 <= 1
 f915 <= 1
 f910 <= 1
 f912 <= 1
 f913 <= 1
 f916 <= 1
 f917 <= 1
 f919 <= 1
 f341 <= 1
 f344 <= 1
 f345 <= 1
 f340 <= 1
 f342 <= 1
 f343 <= 1
 f346 <= 1
 f347 <= 1
 f349 <= 1
 f431 <= 1
 f434 <= 1
 f435 <= 1
 f430 <= 1
 f432 <= 1
 f433 <= 1
 f436 <= 1
 f437 <= 1
 f439 <= 1
 f351 <= 1
 f354 <= 1
 f355 <= 1
 f350 <= 1
 f352 <= 1
 f353 <= 1
 f356 <= 1
 f357 <= 1
 f359 <= 1
 f531 <= 1
 f534 <= 1
 f535 <= 1
 f530 <= 1
 f532 <= 1
 f533 <= 1
 f536 <= 1
 f537 <= 1
 f539 <= 1
 f371 <= 1
 f374 <= 1
 f375 <= 1
 f370 <= 1
 f372 <= 1
 f373 <= 1
 f376 <= 1
 f377 <= 1
 f379 <= 1
 f731 <= 1
 f734 <= 1
 f735 <= 1
 f730 <= 1
 f732 <= 1
 f733 <= 1
 f736 <= 1
 f737 <= 1
 f739 <= 1
 f381 <= 1
 f384 <= 1
 f385 <= 1
 f380 <= 1
 f382 <= 1
 f383 <= 1
 f386 <= 1
 f387 <= 1
 f389 <= 1
 f831 <= 1
 f834 <= 1
 f835 <= 1
 f830 <= 1
 f832 <= 1
 f833 <= 1
 f836 <= 1
 f837 <= 1
 f839 <= 1
 f391 <= 1
 f394 <= 1
 f395 <= 1
 f390 <= 1
 f392 <= 1
 f393 <= 1
 f396 <= 1
 f397 <= 1
 f399 <= 1
 f931 <= 1
 f934 <= 1
 f935 <= 1
 f930 <= 1
 f932 <= 1
 f933 <= 1
 f936 <= 1
 f937 <= 1
 f939 <= 1
 f451 <= 1
 f454 <= 1
 f455 <= 1
 f450 <= 1
 f452 <= 1
 f453 <= 1
 f456 <= 1
 f457 <= 1
 f459 <= 1
 f541 <= 1
 f544 <= 1
 f545 <= 1
 f540 <= 1
 f542 <= 1
 f543 <= 1
 f546 <= 1
 f547 <= 1
 f549 <= 1
 f471 <= 1
 f474 <= 1
 f475 <= 1
 f470 <= 1
 f472 <= 1
 f473 <= 1
 f476 <= 1
 f477 <= 1
 f479 <= 1
 f741 <= 1
 f744 <= 1
 f745 <= 1
 f740 <= 1
 f742 <= 1
 f743 <= 1
 f746 <= 1
 f747 <= 1
 f749 <= 1
 f491 <= 1
 f494 <= 1
 f495 <= 1
 f490 <= 1
 f492 <= 1
 f493 <= 1
 f496 <= 1
 f497 <= 1
 f499 <= 1
 f941 <= 1
 f944 <= 1
 f945 <= 1
 f940 <= 1
 f942 <= 1
 f943 <= 1
 f946 <= 1
 f947 <= 1
 f949 <= 1
 f561 <= 1
 f564 <= 1
 f565 <= 1
 f560 <= 1
 f562 <= 1
 f563 <= 1
 f566 <= 1
 f567 <= 1
 f569 <= 1
 f651 <= 1
 f654 <= 1
 f655 <= 1
 f650 <= 1
 f652 <= 1
 f653 <= 1
 f656 <= 1
 f657 <= 1
 f659 <= 1
 f571 <= 1
 f574 <= 1
 f575 <= 1
 f570 <= 1
 f572 <= 1
 f573 <= 1
 f576 <= 1
 f577 <= 1
 f579 <= 1
 f751 <= 1
 f754 <= 1
 f755 <= 1
 f750 <= 1
 f752 <= 1
 f753 <= 1
 f756 <= 1
 f757 <= 1
 f759 <= 1
 f581 <= 1
 f584 <= 1
 f585 <= 1
 f580 <= 1
 f582 <= 1
 f583 <= 1
 f586 <= 1
 f587 <= 1
 f589 <= 1
 f851 <= 1
 f854 <= 1
 f855 <= 1
 f850 <= 1
 f852 <= 1
 f853 <= 1
 f856 <= 1
 f857 <= 1
 f859 <= 1
 f591 <= 1
 f594 <= 1
 f595 <= 1
 f590 <= 1
 f592 <= 1
 f593 <= 1
 f596 <= 1
 f597 <= 1
 f599 <= 1
 f951 <= 1
 f954 <= 1
 f955 <= 1
 f950 <= 1
 f952 <= 1
 f953 <= 1
 f956 <= 1
 f957 <= 1
 f959 <= 1
 f681 <= 1
 f684 <= 1
 f685 <= 1
 f680 <= 1
 f682 <= 1
 f683 <= 1
 f686 <= 1
 f687 <= 1
 f689 <= 1
 f861 <= 1
 f864 <= 1
 f865 <= 1
 f860 <= 1
 f862 <= 1
 f863 <= 1
 f866 <= 1
 f867 <= 1
 f869 <= 1
 f791 <= 1
 f794 <= 1
 f795 <= 1
 f790 <= 1
 f792 <= 1
 f793 <= 1
 f796 <= 1
 f797 <= 1
 f799 <= 1
 f971 <= 1
 f974 <= 1
 f975 <= 1
 f970 <= 1
 f972 <= 1
 f973 <= 1
 f976 <= 1
 f977 <= 1
 f979 <= 1
 f891 <= 1
 f894 <= 1
 f895 <= 1
 f890 <= 1
 f892 <= 1
 f893 <= 1
 f896 <= 1
 f897 <= 1
 f899 <= 1
 f981 <= 1
 f984 <= 1
 f985 <= 1
 f980 <= 1
 f982 <= 1
 f983 <= 1
 f986 <= 1
 f987 <= 1
 f989 <= 1
Binaries
 y03 y30 y05 y50 y06 y60 y08 y80 y14 y41 y15 y51 y17 y71 y19 y91 y34 y43
 y35 y53 y37 y73 y38 y83 y39 y93 y45 y54 y47 y74 y49 y94 y56 y65 y57 y75
 y58 y85 y59 y95 y68 y86 y79 y97 y89 y98 z0 z1 z2 z3 z4 z5 z6 z7 z8 z9
End
