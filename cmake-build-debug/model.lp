\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  z0 + 0 z1 + 0 z2 + 0 z3 + 0 z4 + 0 z5 + 0 z6 + z7 + 0 z8 + 0 z9
Subject To
 root_flow: - f030 - f130 + f300 + f310 + f350 + f390 - f530 - f930 = 1
 root_flow: - f037 - f137 + f307 + f317 + f357 + f397 - f537 - f937 = 1
 flow_conservation1: f010 - f100 - f130 - f190 + f310 + f910 = 0
 flow_conservation2: = 0
 flow_conservation4: = 0
 flow_conservation5: f350 - f530 - f570 + f750 = 0
 flow_conservation6: = 0
 flow_conservation7: f570 - f750 = 0
 flow_conservation8: = 0
 flow_conservation9: f090 + f190 + f390 - f900 - f910 - f930 = 0
 flow_conservation0: - f017 - f037 - f097 + f107 + f307 + f907 = 0
 flow_conservation1: f017 - f107 - f137 - f197 + f317 + f917 = 0
 flow_conservation2: = 0
 flow_conservation4: = 0
 flow_conservation5: f357 - f537 - f577 + f757 = 0
 flow_conservation6: = 0
 flow_conservation8: = 0
 flow_conservation9: f097 + f197 + f397 - f907 - f917 - f937 = 0
 flow_conservation0: - f011 - f031 - f091 + f101 + f301 + f901 = 0
 flow_conservation2: = 0
 flow_conservation4: = 0
 flow_conservation5: f351 - f531 - f571 + f751 = 0
 flow_conservation6: = 0
 flow_conservation7: f571 - f751 = 0
 flow_conservation8: = 0
 flow_conservation9: f091 + f191 + f391 - f901 - f911 - f931 = 0
 flow_conservation0: - f015 - f035 - f095 + f105 + f305 + f905 = 0
 flow_conservation1: f015 - f105 - f135 - f195 + f315 + f915 = 0
 flow_conservation2: = 0
 flow_conservation4: = 0
 flow_conservation6: = 0
 flow_conservation7: f575 - f755 = 0
 flow_conservation8: = 0
 flow_conservation9: f095 + f195 + f395 - f905 - f915 - f935 = 0
 flow_conservation0: - f019 - f039 - f099 + f109 + f309 + f909 = 0
 flow_conservation1: f019 - f109 - f139 - f199 + f319 + f919 = 0
 flow_conservation2: = 0
 flow_conservation4: = 0
 flow_conservation5: f359 - f539 - f579 + f759 = 0
 flow_conservation6: = 0
 flow_conservation7: f579 - f759 = 0
 flow_conservation8: = 0
 flow_on_terminals: f010 + f030 + f090 - f100 - f300 - f900 = -1
 flow_on_terminals: - f577 + f757 = -1
 flow_on_terminals: - f011 + f101 + f131 + f191 - f311 - f911 = -1
 flow_on_terminals: - f355 + f535 + f575 - f755 = -1
 flow_on_terminals: - f099 - f199 - f399 + f909 + f919 + f939 = -1
 f_and_y_relation: - y01 + f010 <= 0
 f_and_y_relation: - y01 + f017 <= 0
 f_and_y_relation: - y01 + f011 <= 0
 f_and_y_relation: - y01 + f015 <= 0
 f_and_y_relation: - y01 + f019 <= 0
 f_and_y_relation: - y03 + f030 <= 0
 f_and_y_relation: - y03 + f037 <= 0
 f_and_y_relation: - y03 + f031 <= 0
 f_and_y_relation: - y03 + f035 <= 0
 f_and_y_relation: - y03 + f039 <= 0
 f_and_y_relation: - y09 + f090 <= 0
 f_and_y_relation: - y09 + f097 <= 0
 f_and_y_relation: - y09 + f091 <= 0
 f_and_y_relation: - y09 + f095 <= 0
 f_and_y_relation: - y09 + f099 <= 0
 f_and_y_relation: - y10 + f100 <= 0
 f_and_y_relation: - y10 + f107 <= 0
 f_and_y_relation: - y10 + f101 <= 0
 f_and_y_relation: - y10 + f105 <= 0
 f_and_y_relation: - y10 + f109 <= 0
 f_and_y_relation: - y13 + f130 <= 0
 f_and_y_relation: - y13 + f137 <= 0
 f_and_y_relation: - y13 + f131 <= 0
 f_and_y_relation: - y13 + f135 <= 0
 f_and_y_relation: - y13 + f139 <= 0
 f_and_y_relation: - y19 + f190 <= 0
 f_and_y_relation: - y19 + f197 <= 0
 f_and_y_relation: - y19 + f191 <= 0
 f_and_y_relation: - y19 + f195 <= 0
 f_and_y_relation: - y19 + f199 <= 0
 f_and_y_relation: - y30 + f300 <= 0
 f_and_y_relation: - y30 + f307 <= 0
 f_and_y_relation: - y30 + f301 <= 0
 f_and_y_relation: - y30 + f305 <= 0
 f_and_y_relation: - y30 + f309 <= 0
 f_and_y_relation: - y31 + f310 <= 0
 f_and_y_relation: - y31 + f317 <= 0
 f_and_y_relation: - y31 + f311 <= 0
 f_and_y_relation: - y31 + f315 <= 0
 f_and_y_relation: - y31 + f319 <= 0
 f_and_y_relation: - y35 + f350 <= 0
 f_and_y_relation: - y35 + f357 <= 0
 f_and_y_relation: - y35 + f351 <= 0
 f_and_y_relation: - y35 + f355 <= 0
 f_and_y_relation: - y35 + f359 <= 0
 f_and_y_relation: - y39 + f390 <= 0
 f_and_y_relation: - y39 + f397 <= 0
 f_and_y_relation: - y39 + f391 <= 0
 f_and_y_relation: - y39 + f395 <= 0
 f_and_y_relation: - y39 + f399 <= 0
 f_and_y_relation: - y53 + f530 <= 0
 f_and_y_relation: - y53 + f537 <= 0
 f_and_y_relation: - y53 + f531 <= 0
 f_and_y_relation: - y53 + f535 <= 0
 f_and_y_relation: - y53 + f539 <= 0
 f_and_y_relation: - y57 + f570 <= 0
 f_and_y_relation: - y57 + f577 <= 0
 f_and_y_relation: - y57 + f571 <= 0
 f_and_y_relation: - y57 + f575 <= 0
 f_and_y_relation: - y57 + f579 <= 0
 f_and_y_relation: - y75 + f750 <= 0
 f_and_y_relation: - y75 + f757 <= 0
 f_and_y_relation: - y75 + f751 <= 0
 f_and_y_relation: - y75 + f755 <= 0
 f_and_y_relation: - y75 + f759 <= 0
 f_and_y_relation: - y90 + f900 <= 0
 f_and_y_relation: - y90 + f907 <= 0
 f_and_y_relation: - y90 + f901 <= 0
 f_and_y_relation: - y90 + f905 <= 0
 f_and_y_relation: - y90 + f909 <= 0
 f_and_y_relation: - y91 + f910 <= 0
 f_and_y_relation: - y91 + f917 <= 0
 f_and_y_relation: - y91 + f911 <= 0
 f_and_y_relation: - y91 + f915 <= 0
 f_and_y_relation: - y91 + f919 <= 0
 f_and_y_relation: - y93 + f930 <= 0
 f_and_y_relation: - y93 + f937 <= 0
 f_and_y_relation: - y93 + f931 <= 0
 f_and_y_relation: - y93 + f935 <= 0
 f_and_y_relation: - y93 + f939 <= 0
 maximum_of_arcs: y01 + y03 + y09 + y10 + y13 + y19 + y30 + y31 + y35 + y39
   + y53 + y57 + y75 + y90 + y91 + y93 = 9
 delay_limit: 2098 f010 + 2109 f030 + 2109 f090 + 2098 f100 + 1962 f130
   + 2052 f190 + 2109 f300 + 1962 f310 + 1829 f350 + 1950 f390 + 1829 f530
   + 9745 f570 + 9745 f750 + 2109 f900 + 2052 f910 + 1950 f930 - 16158 z0
   <= 3624
 jitter_limit: 4264 f010 + 4265 f030 + 4265 f090 + 4264 f100 + 3562 f130
   + 339 f190 + 4265 f300 + 3562 f310 + 817 f350 + 1102 f390 + 817 f530
   + 2239 f570 + 2239 f750 + 4265 f900 + 339 f910 + 1102 f930 - 11353 z0
   <= 4765
 delay_limit: 2098 f017 + 2109 f037 + 2109 f097 + 2098 f107 + 1962 f137
   + 2052 f197 + 2109 f307 + 1962 f317 + 1829 f357 + 1950 f397 + 1829 f537
   + 9745 f577 + 9745 f757 + 2109 f907 + 2052 f917 + 1950 f937 - 16158 z7
   <= 3624
 jitter_limit: 4264 f017 + 4265 f037 + 4265 f097 + 4264 f107 + 3562 f137
   + 339 f197 + 4265 f307 + 3562 f317 + 817 f357 + 1102 f397 + 817 f537
   + 2239 f577 + 2239 f757 + 4265 f907 + 339 f917 + 1102 f937 - 11353 z7
   <= 4765
 limit_of_variation_between_pairs: 2098 f010 - 2098 f017 + 2109 f030
   - 2109 f037 + 2109 f090 - 2109 f097 + 2098 f100 - 2098 f107 + 1962 f130
   - 1962 f137 + 2052 f190 - 2052 f197 + 2109 f300 - 2109 f307 + 1962 f310
   - 1962 f317 + 1829 f350 - 1829 f357 + 1950 f390 - 1950 f397 + 1829 f530
   - 1829 f537 + 9745 f570 - 9745 f577 + 9745 f750 - 9745 f757 + 2109 f900
   - 2109 f907 + 2052 f910 - 2052 f917 + 1950 f930 - 1950 f937 - 16158 z0
   + 9465 z7 <= 1515
 limit_of_variation_between_pairs: - 2098 f010 + 2098 f017 - 2109 f030
   + 2109 f037 - 2109 f090 + 2109 f097 - 2098 f100 + 2098 f107 - 1962 f130
   + 1962 f137 - 2052 f190 + 2052 f197 - 2109 f300 + 2109 f307 - 1962 f310
   + 1962 f317 - 1829 f350 + 1829 f357 - 1950 f390 + 1950 f397 - 1829 f530
   + 1829 f537 - 9745 f570 + 9745 f577 - 9745 f750 + 9745 f757 - 2109 f900
   + 2109 f907 - 2052 f910 + 2052 f917 - 1950 f930 + 1950 f937 - 16158 z7
   <= 1515
Bounds
Binaries
 y01 f010 f017 f011 f015 f019 y03 f030 f037 f031 f035 f039 y09 f090 f097
 f091 f095 f099 y10 f100 f107 f101 f105 f109 y13 f130 f137 f131 f135 f139
 y19 f190 f197 f191 f195 f199 y30 f300 f307 f301 f305 f309 y31 f310 f317
 f311 f315 f319 y35 f350 f357 f351 f355 f359 y39 f390 f397 f391 f395 f399
 y53 f530 f537 f531 f535 f539 y57 f570 f577 f571 f575 f579 y75 f750 f757
 f751 f755 f759 y90 f900 f907 f901 f905 f909 y91 f910 f917 f911 f915 f919
 y93 f930 f937 f931 f935 f939 z0 z1 z2 z3 z4 z5 z6 z7 z8 z9
End
