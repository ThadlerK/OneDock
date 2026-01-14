from pymol import cmd,stored

set depth_cue, 1
set fog_start, 0.4

set_color b_col, [36,36,85]
set_color t_col, [10,10,10]
set bg_rgb_bottom, b_col
set bg_rgb_top, t_col      
set bg_gradient

set  spec_power  =  200
set  spec_refl   =  0

load "data/8WM3_Tiagabine_ACE2.pdb", protein
create ligands, protein and organic
select xlig, protein and organic
delete xlig

hide everything, all

color white, elem c
color bluewhite, protein
#show_as cartoon, protein
show surface, protein
#set transparency, 0.15

show sticks, ligands
set stick_color, magenta




# SAS points

load "data/8WM3_Tiagabine_ACE2.pdb_points.pdb.gz", points
hide nonbonded, points
show nb_spheres, points
set sphere_scale, 0.2, points
cmd.spectrum("b", "green_red", selection="points", minimum=0, maximum=0.7)


stored.list=[]
cmd.iterate("(resn STP)","stored.list.append(resi)")    # read info about residues STP
lastSTP=stored.list[-1] # get the index of the last residue
hide lines, resn STP

cmd.select("rest", "resn STP and resi 0")

for my_index in range(1,int(lastSTP)+1): cmd.select("pocket"+str(my_index), "resn STP and resi "+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.show("spheres","pocket"+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.set("sphere_scale","0.4","pocket"+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.set("sphere_transparency","0.1","pocket"+str(my_index))



set_color pcol1 = [0.361,0.576,0.902]
select surf_pocket1, protein and id [17897,17323,17307,17888,17893,17896,17934,17920,17887,18259,18261,18256,18301,18302,18303,17393,17349,17350,18338,17392,17925,17461,17463,18324,18325,17923,17460,17304] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.278,0.412,0.702]
select surf_pocket2, protein and id [10664,10666,10668,10667,10643,10644,12579,12580,12611,12608,12609,12610,10671,12710,10674,12711,12681,12680,12682,12683,12641,12638,10617,10621,10622,10619,10673,12554,12556,11001,12574,12575,12577,12501,12507,12555,12576,11000,12511,12645,12505] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.361,0.486,0.902]
select surf_pocket3, protein and id [10989,10990,14232,14233,14235,14236,11012,11092,12520,12524,12525,12522,14212,12527,14204,11104,11105,11179,12539,12540,10992,12563,12561,12564,11203,11204,11074,11075,11071,10955] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.278,0.345,0.702]
select surf_pocket4, protein and id [14514,10983,14515,14262,14312,14311,14398,14349,14350,14351,14433,14434,14424,14439,14443,14471,14483,11006,11026,14392,14395,14396,10973,10974,10975,12654,12637,14394,14404,14413,14393,11003,11004,11005,11007,11031,11013,14353] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.361,0.400,0.902]
select surf_pocket5, protein and id [17987,18202,18015,18225,18226,18229,18220,18232,19470,19471,19436,19365,19251,17540,17528,18010,17507,18224,18219] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.282,0.278,0.702]
select surf_pocket6, protein and id [11012,11092,14344,14345,11132,11013,14353,14235,14236,11105,11130,11131,10983,14263,14515,14258,14262,14302,14312,14301,14350,14351,14292,14342,14291,14254,14256,14257,14255,14293] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.408,0.361,0.902]
select surf_pocket7, protein and id [10694,10683,12499,12501,12506,12507,12512,12555,12554,10686,13709,13703,10688,10695,12496,12497,12490,11297,12495,12537,13737,13705,13708,11298,13707,13735,13736,12513,12515,12535,12547,12548,13727,10673,10678] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.353,0.278,0.702]
select surf_pocket8, protein and id [14705,14706,14707,14658,14659,14660,14661,14662,14663,14664,14655,12404,12405,12402,12433,14626,14628,14630,14625,14629,14627,14601,14623,14632,14633,14636,14648,14649,14654,14650,14622,14138,14140,14141,14143,14147,14126,14132,14166,14884,14904,14903,14932,12440] 
set surface_color,  pcol8, surf_pocket8 
set_color pcol9 = [0.498,0.361,0.902]
select surf_pocket9, protein and id [18852,18853,18854,18907,20431,19879,20430,19909,19910,19926,20423,19928,18926,19884,19886,18775,18778,18779,18780,20424,20420,20415,18781,18825,18828,18810,18774,18777,18812,18932,18940,18929,18941,18947] 
set surface_color,  pcol9, surf_pocket9 
set_color pcol10 = [0.420,0.278,0.702]
select surf_pocket10, protein and id [14733,14734,14729,10797,10798,14711,14737,14718,14721,13420,13421,13424,13429,14732,14095,14093,14094,14124,14125,14126,14127,14130,14132,14135,14138,10809,10775,12441,12440,14092,14091,14129,14090,12404,12402,12433,10807,13415,13416,13418] 
set surface_color,  pcol10, surf_pocket10 
set_color pcol11 = [0.584,0.361,0.902]
select surf_pocket11, protein and id [11397,11398,11399,13670,11406,11408,11361,11363,11365,21000,21030,21029,13697,11943,11944,21057,21054,13604,13643,13645,11436,11437] 
set surface_color,  pcol11, surf_pocket11 
set_color pcol12 = [0.490,0.278,0.702]
select surf_pocket12, protein and id [16618,19180,19537,19538,16765,16700,16701,16739,19631] 
set surface_color,  pcol12, surf_pocket12 
set_color pcol13 = [0.675,0.361,0.902]
select surf_pocket13, protein and id [11074,11070,11081,11083,11067,11076,11080,11204,12531,11211,13775,13777,13784,13807,13753,12566,12590,11059,11066,12540,12542,13748,12569,13776,12568] 
set surface_color,  pcol13, surf_pocket13 
set_color pcol14 = [0.557,0.278,0.702]
select surf_pocket14, protein and id [13759,13760,13761,13744,13770,13930,11266,11267,11928,13960,11233,13937] 
set surface_color,  pcol14, surf_pocket14 
set_color pcol15 = [0.761,0.361,0.902]
select surf_pocket15, protein and id [11224,14159,11222,11247,11250,11220,11225,14105,14107,14924,14146,14925,14173,14149,14152,14153,14154] 
set surface_color,  pcol15, surf_pocket15 
set_color pcol16 = [0.627,0.278,0.702]
select surf_pocket16, protein and id [18803,18844,18845,16249,18818,16291,18890,17082,17081,18990,18983,18984,18985,18989,18804,18597,19011,16292,17124,17132,17133,17135,18568,18570,19014,18602,18604] 
set surface_color,  pcol16, surf_pocket16 
set_color pcol17 = [0.851,0.361,0.902]
select surf_pocket17, protein and id [14587,14618,14616,14557,14582,14561,12418,12420,12419,12423,12421,10889,12424,10888,14560,12395,14589,10896,10897,14590,10926,10927,12394,12417,12422] 
set surface_color,  pcol17, surf_pocket17 
set_color pcol18 = [0.698,0.278,0.702]
select surf_pocket18, protein and id [17964,17997,17753,17754,17994,17996,17992,17993,18034,17769,17773,18179,18181,18182,18033,17776,18116,17739,17744,17748,17750,17974,17975,17735,18177,18178,18180,18130,18188,18191] 
set surface_color,  pcol18, surf_pocket18 
set_color pcol19 = [0.902,0.361,0.859]
select surf_pocket19, protein and id [20934,11681,11727,20956,11725,11458,11471,11474,11475,11442,11443,11445,20929,11477,11724,20935,11516,11728,11729,11500,11417,11419,11377,20977,11860] 
set surface_color,  pcol19, surf_pocket19 
set_color pcol20 = [0.702,0.278,0.635]
select surf_pocket20, protein and id [19914,19948,20575,17018,17020,19942,19915,17013,18911,20593,18920,20591,20568,16182,16184,16221] 
set surface_color,  pcol20, surf_pocket20 
set_color pcol21 = [0.902,0.361,0.773]
select surf_pocket21, protein and id [18925,18971,18973,20380,20407,20409,20381,18946,16328,16329,20378,18921,16300,16309,20571,16260,16263] 
set surface_color,  pcol21, surf_pocket21 
set_color pcol22 = [0.702,0.278,0.565]
select surf_pocket22, protein and id [14732,14050,14059,14095,14093,14094,13421,11696,13424,13429,14733,14734,14739,14740,14743,14747,14805,11706] 
set surface_color,  pcol22, surf_pocket22 
set_color pcol23 = [0.902,0.361,0.682]
select surf_pocket23, protein and id [11618,11571,11572,11573,11574,11575,11576,11608,11570,11504,11643,11644,13523,11494,11637,11638,11639,11508,11609,11607,11634] 
set surface_color,  pcol23, surf_pocket23 
set_color pcol24 = [0.702,0.278,0.498]
select surf_pocket24, protein and id [17155,17166,17168,19119,19117,18544,18515,18518,18526,19121,19160,19152,18190,18216,18214,19189] 
set surface_color,  pcol24, surf_pocket24 
set_color pcol25 = [0.902,0.361,0.596]
select surf_pocket25, protein and id [11239,11241,11242,13792,13951,13762,13791,13926,13930,13931,11240] 
set surface_color,  pcol25, surf_pocket25 
set_color pcol26 = [0.702,0.278,0.427]
select surf_pocket26, protein and id [14073,14075,14114,14897,14924,13997,13975,13998,11251,11253,13979,13977,14113,11274] 
set surface_color,  pcol26, surf_pocket26 
set_color pcol27 = [0.902,0.361,0.506]
select surf_pocket27, protein and id [10681,12118,10702,12765,12152,13655,10710,10706,10711] 
set surface_color,  pcol27, surf_pocket27 
set_color pcol28 = [0.702,0.278,0.361]
select surf_pocket28, protein and id [10719,10721,11320,11322,12497,11295,11293,11294,11296,11298,10689,10690,10693,14121,11284] 
set surface_color,  pcol28, surf_pocket28 
set_color pcol29 = [0.902,0.361,0.420]
select surf_pocket29, protein and id [14475,14477,14502,12664,14476,10913,10909,10942,10943,10939,10941,12695,12661,10936,10938] 
set surface_color,  pcol29, surf_pocket29 
set_color pcol30 = [0.702,0.278,0.290]
select surf_pocket30, protein and id [15523,15733,15734,18103,15742,18089,18092,18094,15524,15492,18104,18098] 
set surface_color,  pcol30, surf_pocket30 
set_color pcol31 = [0.902,0.388,0.361]
select surf_pocket31, protein and id [14589,12451,10896,14590,10926,10927,14569,14188,14557,12484,12486] 
set surface_color,  pcol31, surf_pocket31 
set_color pcol32 = [0.702,0.337,0.278]
select surf_pocket32, protein and id [19219,19221,18170,19560,19561,19461,19557,19558,19589,19559,19567,19565,19175,19177] 
set surface_color,  pcol32, surf_pocket32 
set_color pcol33 = [0.902,0.478,0.361]
select surf_pocket33, protein and id [20417,20183,20184,20116,20117,18943,20118,18946,18944,18945,20115,20141,20142,20438,20144] 
set surface_color,  pcol33, surf_pocket33 
set_color pcol34 = [0.702,0.408,0.278]
select surf_pocket34, protein and id [12430,12432,14627,14187,14599,14601,14589,12451,12427,14166,14163,14164,14189,12457] 
set surface_color,  pcol34, surf_pocket34 
set_color pcol35 = [0.902,0.565,0.361]
select surf_pocket35, protein and id [19020,19022,15951,19019,16077,19034,17153,17155,17120,17119,17141,16112,19013] 
set surface_color,  pcol35, surf_pocket35 
set_color pcol36 = [0.702,0.475,0.278]
select surf_pocket36, protein and id [17817,17819,18041,18042,18039,17569,17570,17801,17589] 
set surface_color,  pcol36, surf_pocket36 
set_color pcol37 = [0.902,0.655,0.361]
select surf_pocket37, protein and id [18845,18890,18891,18879,18566,17146,18537,17132,17144,17106] 
set surface_color,  pcol37, surf_pocket37 
set_color pcol38 = [0.702,0.545,0.278]
select surf_pocket38, protein and id [12905,12907,12873,12875,12901,12902,12903,12904,13523,13522,13501,13531,12874,12877] 
set surface_color,  pcol38, surf_pocket38 
set_color pcol39 = [0.902,0.741,0.361]
select surf_pocket39, protein and id [10870,10872,12475,12742,12743,12467,12737,12771,12738,12741,10704,12502,12503] 
set surface_color,  pcol39, surf_pocket39 
set_color pcol40 = [0.702,0.612,0.278]
select surf_pocket40, protein and id [18714,18760,18726,18698,18699,18727,18730,18734,18763,16420,16421,18674,16428,16453] 
set surface_color,  pcol40, surf_pocket40 
set_color pcol41 = [0.902,0.831,0.361]
select surf_pocket41, protein and id [18458,18460,18462,17305,17308,18276,17300,17301,17304,17307,18261,18262,18270,18271] 
set surface_color,  pcol41, surf_pocket41 
set_color pcol42 = [0.702,0.682,0.278]
select surf_pocket42, protein and id [18326,17521,17522,18320,18316,18317,18318,19370,19371,19357,18293] 
set surface_color,  pcol42, surf_pocket42 
set_color pcol43 = [0.878,0.902,0.361]
select surf_pocket43, protein and id [10709,10714,10717,13662,13663,10723,11321,10712,13447] 
set surface_color,  pcol43, surf_pocket43 
set_color pcol44 = [0.651,0.702,0.278]
select surf_pocket44, protein and id [19218,19214,19246,19242,19243,19244,19590,19435,19424,19429] 
set surface_color,  pcol44, surf_pocket44 
   

deselect

orient
