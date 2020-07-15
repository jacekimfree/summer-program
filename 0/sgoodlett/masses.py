
charge = {
"X"  :0  ,"H"  :1  ,"HE" :2  ,"LI" :3  ,"BE" :4  ,"B"  :5  ,"C"  :6  ,
"N"  :7  ,"O"  :8  ,"F"  :9  ,"NE" :10 ,"NA" :11 ,"MG" :12 ,"AL" :13 ,
"SI" :14 ,"P"  :15 ,"S"  :16 ,"CL" :17 ,"AR" :18 ,"K"  :19 ,"CA" :20 ,
"SC" :21 ,"TI" :22 ,"V"  :23 ,"CR" :24 ,"MN" :25 ,"FE" :26 ,"CO" :27 ,
"NI" :28 ,"CU" :29 ,"ZN" :30 ,"GA" :31 ,"GE" :32 ,"AS" :33 ,"SE" :34 ,
"BR" :35 ,"KR" :36 ,"RB" :37 ,"SR" :38 ,"Y"  :39 ,"ZR" :40 ,"NB" :41 ,
"MO" :42 ,"TC" :43 ,"RU" :44 ,"RH" :45 ,"PD" :46 ,"AG" :47 ,"CD" :48 ,
"IN" :49 ,"SN" :50 ,"SB" :51 ,"TE" :52 ,"I"  :53 ,"XE" :54 ,"CS" :55 ,
"BA" :56 ,"LA" :57 ,"CE" :58 ,"PR" :59 ,"ND" :60 ,"PM" :61 ,"SM" :62 ,
"EU" :63 ,"GD" :64 ,"TB" :65 ,"DY" :66 ,"HO" :67 ,"ER" :68 ,"TM" :69 ,
"YB" :70 ,"LU" :71 ,"HF" :72 ,"TA" :73 ,"W"  :74 ,"RE" :75 ,"OS" :76 ,
"IR" :77 ,"PT" :78 ,"AU" :79 ,"HG" :80 ,"TL" :81 ,"PB" :82 ,"BI" :83 ,
"PO" :84 ,"AT" :85 ,"RN" :86 ,"FR" :87 ,"RA" :88 ,"AC" :89 ,"TH" :90 ,
"PA" :91 ,"U"  :92 ,"NP" :93 ,"PU" :94 ,"AM" :95 ,"CM" :96 ,"BK" :97 ,
"CF" :98 ,"ES" :99 ,"FM" :100,"MD" :101,"NO" :102,"LR" :103,"RF" :104,
"DB" :105,"SG" :106,"BH" :107,"HS" :108,"MT" :109,"DS" :110,"RG" :111,
"UUB":112,"UUT":113,"UUQ":114,"UUP":115,"UUH":116,"UUS":117,"UUO":118 }

mass = [
0.0,1.00782503207,4.00260325415,7.016004548,9.012182201,11.009305406,
12,14.00307400478,15.99491461956,18.998403224,19.99244017542,
22.98976928087,23.985041699,26.981538627,27.97692653246,30.973761629,
31.972070999,34.968852682,39.96238312251,38.963706679,39.962590983,
44.955911909,47.947946281,50.943959507,51.940507472,54.938045141,
55.934937475,58.933195048,57.935342907,62.929597474,63.929142222,
68.925573587,73.921177767,74.921596478,79.916521271,78.918337087,
85.910610729,84.911789737,87.905612124,88.905848295,89.904704416,
92.906378058,97.905408169,98.906254747,101.904349312,102.905504292,
105.903485715,106.90509682,113.90335854,114.903878484,119.902194676,
120.903815686,129.906224399,126.904472681,131.904153457,132.905451932,
137.905247237,138.906353267,139.905438706,140.907652769,141.907723297,
144.912749023,151.919732425,152.921230339,157.924103912,158.925346757,
163.929174751,164.93032207,165.930293061,168.93421325,173.938862089,
174.940771819,179.946549953,180.947995763,183.950931188,186.955753109,
191.96148069,192.96292643,194.964791134,196.966568662,201.970643011,
204.974427541,207.976652071,208.980398734,208.982430435,210.987496271,
222.017577738,222.01755173,228.031070292,227.027752127,232.038055325,
231.03588399,238.050788247,237.048173444,242.058742611,243.06138108,
247.07035354,247.07030708,251.079586788,252.082978512,257.095104724,
258.098431319,255.093241131,260.105504,263.112547,255.107398,259.114500,
262.122892,263.128558,265.136151,281.162061,272.153615,283.171792,283.176451,
285.183698,287.191186,292.199786,291.206564,293.214670]

def get_mass(atom):
    return mass[charge[atom.upper()]]

def get_charge(atom):
    return charge[atom.upper()]