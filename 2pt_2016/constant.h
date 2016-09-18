#ifndef CONSTANT_H
#define CONSTANT_H
#include <cmath>

template <class SW> SW SWAP(SW *a, SW *b)
{
        SW swap=*a;
        *a=*b;
        *b=swap;
        return 0;
}

const int    CPNAME=128;                //maximum length of a compound
const int    maxcnt=6;                  //maximum number of connectivity per atom
const int    maxffp=3;                  //maximum number of force field parameters
const double MEGABYTE=1024*1024;        // 
const double Na=6.0221367E23;           //Avogadro's number
const double epslon0=8.8542E-12;        //permitivity C2/N/m2
const double e2C=1.602188E-19;          //C per electron
const double eleu2kcalpmol=1.0/(4.0*3.14159265358979323846*epslon0)*Na*(1.602188E-19*1.602188E-19)*1E10/4184; //=332.007 kcal A/e2/mol
const double tolerance=1e-14;           //tolerence
const double dblprecision=1e-16;
const double PI=3.141592653589793116;
const double eA2debye=4.802;            // 1e x 1A = 4.802 Debye
const double kb=1.380658E-23;           // Boltzmann constant J/K
const double h=6.62606896E-34;          // Plank constant Js
const double R=kb*Na;                   //8.314472;    /gas constant J/mol K
const double vlight=2.99792458E8;        // m/s
const double Bohr2A=0.529177249;        //1 Bohr = 0.529177249 A
const double Hartree2kcalmol=627.5095;   // 1 Hartree = 627.5095 kcal/mol
const double caltoj=4.184;               // 1 cal = 4.184 J

#define true          1
#define false         0
#define none          0     //none
#define filenotopen   1001  //file cannot be opened
#define unregkey      1002  //unrecognized keyword
#define c2trj         2001  //cerius2 trj
#define asctrj        2002  //ASCII trj
#define cpmdtrj       2003  //CPMD trj
#define lmptrj        2004  //LAMMPS trj
#define usrtrj        2005  //USER  trj
#define fgtrj         2006  //fix gromacs trj
#define mol2grp       3001  //create group file for each molecule
#define eqmass2grp    3002  //create group file for indentical masses
#define eqmolmass2grp 3003  //create group file for indentical molecules
#define strtbgf       4001  //read structure from bgf      
#define strtlmp       4002  //read structure from lammps data file
#define strtcpmd      4003  //read structure from CPMD in file
#define strtgro       4004  //read structure from gromacs data file
#define strtg96       4005  //read structure from gromacs G96 data file
#define fflmp         4011  //read focefield parameters from lammps data file
#define ffdreid       4012  //read focefield parameters from Dreiding parameter file
#define EWALD         5001  //Ewald sum
#define SPLINE        5002  //spline
#define DIRECT        5003  //cutoff
#define GEOMETRIC     5004  //geometric mean, sqrt(i*j) for both Do and Ro
#define ARITHMETIC    5005  //arithmetic mean, sqrt(i*j) for Do, (i+j)/2 for Ro
#define SIXTHPOWER    5006  //6th power,
#define TORALL        -001  //analyze each individual torsions
#define TORTYPE       -002  //analyze each torsional type 
#define CLRATOMTYPE      0  //color by atom type
#define CLRFFYPE         1  //color by force field type
#define CLRCHARGE        2  //color by charge
#define CLRVELOCITY      3  //color by velocity (atomic temperature)
#define CLRHIDE          4  //hide everything
#define DRAWID           1  //label atom id
#define DRAWCHARGE       2  //label atom charge
#define DRAWNAME         3  //label atom name

#define DSPLNONE         0  //no display
#define DSPBALL          1  //ball
#define DSPSTICK         2  //stick
#define DSPCYLDR         3  //cylinder
#define DSPSPHERE        4  //sphere

#define VELTYPE                 4
#define vtotal                  VELTYPE
#define vtrans                  VELTYPE-4
#define vangul                  VELTYPE-3
#define vimvib                  VELTYPE-2
#define vrotat                  VELTYPE-1

#define FORCETYPE 9
#define NET FORCETYPE-9
#define BND FORCETYPE-8
#define ANG FORCETYPE-7
#define TOR FORCETYPE-6
#define INV FORCETYPE-5
#define VDW FORCETYPE-4
#define VDL FORCETYPE-3
#define ELS FORCETYPE-2
#define ELL FORCETYPE-1

#define MAPORIGINAL 0
#define MAPDEFAULT  1
#define MAPONECELL  2

#ifdef _WIN32
#define M_PI 3.1415926535897932384626433832795
#define M_SQRT2 1.4142135623730950488016887242097
#endif

#ifdef _WIN32
// On linux, I define this in a Makefile.
// On Windows, I like it defined by default.
#define DEBUG
#endif

#ifdef DEBUG
   #define ASSERT(exp) { if (!(exp)) { \
      printf( "Assertion error at line %d in %s\n", __LINE__, __FILE__ ); \
      /*sprintf(0, "This should cause a crash blah blah blah blah.");*/ \
      exit( 1 ); \
   }}
   #define ASSERT_IS_EQUAL(a,b) { \
      ASSERT( (b) - 0.0005 < (a) && (a) < (b) + 0.0005 ); \
   }
#else
   #define ASSERT(exp)
   #define ASSERT_IS_EQUAL(a,b)
#endif


inline int ROUND( float x ) {
   return x < 0 ? -(int)(-x+0.5f) : (int)(x+0.5f);
}
inline int ROUND( double x ) {
   return x < 0 ? -(int)(-x+0.5f) : (int)(x+0.5f);
}

//used for visualization
#define	HIDE	-1
#define	SNOW	0
#define	GHOSTWHITE	1
#define	WHITESMOKE	2
#define	GAINSBORO	3
#define	FLORALWHITE	4
#define	OLDLACE	5
#define	LINEN	6
#define	ANTIQUEWHITE	7
#define	PAPAYAWHIP	8
#define	BLANCHEDALMOND	9
#define	BISQUE	10
#define	PEACHPUFF	11
#define	NAVAJOWHITE	12
#define	MOCCASIN	13
#define	CORNSILK	14
#define	IVORY	15
#define	LEMONCHIFFON	16
#define	SEASHELL	17
#define	HONEYDEW	18
#define	MINTCREAM	19
#define	AZURE	20
#define	ALICEBLUE	21
#define	LAVENDER	22
#define	LAVENDERBLUSH	23
#define	MISTYROSE	24
#define	WHITE	25
#define	BLACK	26
#define	DARKSLATEGREY	27
#define	DIMGRAY	28
#define	SLATEGREY	29
#define	LIGHTSLATEGREY	30
#define	GREY	31
#define	LIGHTGRAY	32
#define	MIDNIGHTBLUE	33
#define	NAVYBLUE	34
#define	CORNFLOWERBLUE	35
#define	DARKSLATEBLUE	36
#define	SLATEBLUE	37
#define	MEDIUMSLATEBLUE	38
#define	LIGHTSLATEBLUE	39
#define	MEDIUMBLUE	40
#define	ROYALBLUE	41
#define	BLUE	42
#define	DODGERBLUE	43
#define	DEEPSKYBLUE	44
#define	SKYBLUE	45
#define	LIGHTSKYBLUE	46
#define	STEELBLUE	47
#define	LIGHTSTEELBLUE	48
#define	LIGHTBLUE	49
#define	POWDERBLUE	50
#define	PALETURQUOISE	51
#define	DARKTURQUOISE	52
#define	MEDIUMTURQUOISE	53
#define	TURQUOISE	54
#define	CYAN	55
#define	LIGHTCYAN	56
#define	CADETBLUE	57
#define	MEDIUMAQUAMARINE	58
#define	AQUAMARINE	59
#define	DARKGREEN	60
#define	DARKOLIVEGREEN	61
#define	DARKSEAGREEN	62
#define	SEAGREEN	63
#define	MEDIUMSEAGREEN	64
#define	LIGHTSEAGREEN	65
#define	PALEGREEN	66
#define	SPRINGGREEN	67
#define	LAWNGREEN	68
#define	GREEN	69
#define	CHARTREUSE	70
#define	MEDIUMSPRINGGREEN	71
#define	GREENYELLOW	72
#define	LIMEGREEN	73
#define	YELLOWGREEN	74
#define	FORESTGREEN	75
#define	OLIVEDRAB	76
#define	DARKKHAKI	77
#define	KHAKI	78
#define	PALEGOLDENROD	79
#define	LIGHTGOLDENRODYELLOW	80
#define	LIGHTYELLOW	81
#define	YELLOW	82
#define	GOLD	83
#define	LIGHTGOLDENROD	84
#define	GOLDENROD	85
#define	DARKGOLDENROD	86
#define	ROSYBROWN	87
#define	INDIANRED	88
#define	SADDLEBROWN	89
#define	SIENNA	90
#define	PERU	91
#define	BURLYWOOD	92
#define	BEIGE	93
#define	WHEAT	94
#define	SANDYBROWN	95
#define	TAN	96
#define	CHOCOLATE	97
#define	FIREBRICK	98
#define	BROWN	99
#define	DARKSALMON	100
#define	SALMON	101
#define	LIGHTSALMON	102
#define	ORANGE	103
#define	DARKORANGE	104
#define	CORAL	105
#define	LIGHTCORAL	106
#define	TOMATO	107
#define	ORANGERED	108
#define	RED	109
#define	HOTPINK	110
#define	DEEPPINK	111
#define	PINK	112
#define	LIGHTPINK	113
#define	PALEVIOLETRED	114
#define	MAROON	115
#define	MEDIUMVIOLETRED	116
#define	VIOLETRED	117
#define	MAGENTA	118
#define	VIOLET	119
#define	PLUM	120
#define	ORCHID	121
#define	MEDIUMORCHID	122
#define	DARKORCHID	123
#define	DARKVIOLET	124
#define	BLUEVIOLET	125
#define	PURPLE	126
#define	MEDIUMPURPLE	127
#define	THISTLE	128
#define	SNOW1	129
#define	SNOW2	130
#define	SNOW3	131
#define	SNOW4	132
#define	SEASHELL1	133
#define	SEASHELL2	134
#define	SEASHELL3	135
#define	SEASHELL4	136
#define	ANTIQUEWHITE1	137
#define	ANTIQUEWHITE2	138
#define	ANTIQUEWHITE3	139
#define	ANTIQUEWHITE4	140
#define	BISQUE1	141
#define	BISQUE2	142
#define	BISQUE3	143
#define	BISQUE4	144
#define	PEACHPUFF1	145
#define	PEACHPUFF2	146
#define	PEACHPUFF3	147
#define	PEACHPUFF4	148
#define	NAVAJOWHITE1	149
#define	NAVAJOWHITE2	150
#define	NAVAJOWHITE3	151
#define	NAVAJOWHITE4	152
#define	LEMONCHIFFON1	153
#define	LEMONCHIFFON2	154
#define	LEMONCHIFFON3	155
#define	LEMONCHIFFON4	156
#define	CORNSILK1	157
#define	CORNSILK2	158
#define	CORNSILK3	159
#define	CORNSILK4	160
#define	IVORY1	161
#define	IVORY2	162
#define	IVORY3	163
#define	IVORY4	164
#define	HONEYDEW1	165
#define	HONEYDEW2	166
#define	HONEYDEW3	167
#define	HONEYDEW4	168
#define	LAVENDERBLUSH1	169
#define	LAVENDERBLUSH2	170
#define	LAVENDERBLUSH3	171
#define	LAVENDERBLUSH4	172
#define	MISTYROSE1	173
#define	MISTYROSE2	174
#define	MISTYROSE3	175
#define	MISTYROSE4	176
#define	AZURE1	177
#define	AZURE2	178
#define	AZURE3	179
#define	AZURE4	180
#define	SLATEBLUE1	181
#define	SLATEBLUE2	182
#define	SLATEBLUE3	183
#define	SLATEBLUE4	184
#define	ROYALBLUE1	185
#define	ROYALBLUE2	186
#define	ROYALBLUE3	187
#define	ROYALBLUE4	188
#define	BLUE1	189
#define	BLUE2	190
#define	BLUE3	191
#define	BLUE4	192
#define	DODGERBLUE1	193
#define	DODGERBLUE2	194
#define	DODGERBLUE3	195
#define	DODGERBLUE4	196
#define	STEELBLUE1	197
#define	STEELBLUE2	198
#define	STEELBLUE3	199
#define	STEELBLUE4	200
#define	DEEPSKYBLUE1	201
#define	DEEPSKYBLUE2	202
#define	DEEPSKYBLUE3	203
#define	DEEPSKYBLUE4	204
#define	SKYBLUE1	205
#define	SKYBLUE2	206
#define	SKYBLUE3	207
#define	SKYBLUE4	208
#define	LIGHTSKYBLUE1	209
#define	LIGHTSKYBLUE2	210
#define	LIGHTSKYBLUE3	211
#define	LIGHTSKYBLUE4	212
#define	SLATEGRAY1	213
#define	SLATEGRAY2	214
#define	SLATEGRAY3	215
#define	SLATEGRAY4	216
#define	LIGHTSTEELBLUE1	217
#define	LIGHTSTEELBLUE2	218
#define	LIGHTSTEELBLUE3	219
#define	LIGHTSTEELBLUE4	220
#define	LIGHTBLUE1	221
#define	LIGHTBLUE2	222
#define	LIGHTBLUE3	223
#define	LIGHTBLUE4	224
#define	LIGHTCYAN1	225
#define	LIGHTCYAN2	226
#define	LIGHTCYAN3	227
#define	LIGHTCYAN4	228
#define	PALETURQUOISE1	229
#define	PALETURQUOISE2	230
#define	PALETURQUOISE3	231
#define	PALETURQUOISE4	232
#define	CADETBLUE1	233
#define	CADETBLUE2	234
#define	CADETBLUE3	235
#define	CADETBLUE4	236
#define	TURQUOISE1	237
#define	TURQUOISE2	238
#define	TURQUOISE3	239
#define	TURQUOISE4	240
#define	CYAN1	241
#define	CYAN2	242
#define	CYAN3	243
#define	CYAN4	244
#define	DARKSLATEGRAY1	245
#define	DARKSLATEGRAY2	246
#define	DARKSLATEGRAY3	247
#define	DARKSLATEGRAY4	248
#define	AQUAMARINE1	249
#define	AQUAMARINE2	250
#define	AQUAMARINE3	251
#define	AQUAMARINE4	252
#define	DARKSEAGREEN1	253
#define	DARKSEAGREEN2	254
#define	DARKSEAGREEN3	255
#define	DARKSEAGREEN4	256
#define	SEAGREEN1	257
#define	SEAGREEN2	258
#define	SEAGREEN3	259
#define	SEAGREEN4	260
#define	PALEGREEN1	261
#define	PALEGREEN2	262
#define	PALEGREEN3	263
#define	PALEGREEN4	264
#define	SPRINGGREEN1	265
#define	SPRINGGREEN2	266
#define	SPRINGGREEN3	267
#define	SPRINGGREEN4	268
#define	GREEN1	269
#define	GREEN2	270
#define	GREEN3	271
#define	GREEN4	272
#define	CHARTREUSE1	273
#define	CHARTREUSE2	274
#define	CHARTREUSE3	275
#define	CHARTREUSE4	276
#define	OLIVEDRAB1	277
#define	OLIVEDRAB2	278
#define	OLIVEDRAB3	279
#define	OLIVEDRAB4	280
#define	DARKOLIVEGREEN1	281
#define	DARKOLIVEGREEN2	282
#define	DARKOLIVEGREEN3	283
#define	DARKOLIVEGREEN4	284
#define	KHAKI1	285
#define	KHAKI2	286
#define	KHAKI3	287
#define	KHAKI4	288
#define	LIGHTGOLDENROD1	289
#define	LIGHTGOLDENROD2	290
#define	LIGHTGOLDENROD3	291
#define	LIGHTGOLDENROD4	292
#define	LIGHTYELLOW1	293
#define	LIGHTYELLOW2	294
#define	LIGHTYELLOW3	295
#define	LIGHTYELLOW4	296
#define	YELLOW1	297
#define	YELLOW2	298
#define	YELLOW3	299
#define	YELLOW4	300
#define	GOLD1	301
#define	GOLD2	302
#define	GOLD3	303
#define	GOLD4	304
#define	GOLDENROD1	305
#define	GOLDENROD2	306
#define	GOLDENROD3	307
#define	GOLDENROD4	308
#define	DARKGOLDENROD1	309
#define	DARKGOLDENROD2	310
#define	DARKGOLDENROD3	311
#define	DARKGOLDENROD4	312
#define	ROSYBROWN1	313
#define	ROSYBROWN2	314
#define	ROSYBROWN3	315
#define	ROSYBROWN4	316
#define	INDIANRED1	317
#define	INDIANRED2	318
#define	INDIANRED3	319
#define	INDIANRED4	320
#define	SIENNA1	321
#define	SIENNA2	322
#define	SIENNA3	323
#define	SIENNA4	324
#define	BURLYWOOD1	325
#define	BURLYWOOD2	326
#define	BURLYWOOD3	327
#define	BURLYWOOD4	328
#define	WHEAT1	329
#define	WHEAT2	330
#define	WHEAT3	331
#define	WHEAT4	332
#define	TAN1	333
#define	TAN2	334
#define	TAN3	335
#define	TAN4	336
#define	CHOCOLATE1	337
#define	CHOCOLATE2	338
#define	CHOCOLATE3	339
#define	CHOCOLATE4	340
#define	FIREBRICK1	341
#define	FIREBRICK2	342
#define	FIREBRICK3	343
#define	FIREBRICK4	344
#define	BROWN1	345
#define	BROWN2	346
#define	BROWN3	347
#define	BROWN4	348
#define	SALMON1	349
#define	SALMON2	350
#define	SALMON3	351
#define	SALMON4	352
#define	LIGHTSALMON1	353
#define	LIGHTSALMON2	354
#define	LIGHTSALMON3	355
#define	LIGHTSALMON4	356
#define	ORANGE1	357
#define	ORANGE2	358
#define	ORANGE3	359
#define	ORANGE4	360
#define	DARKORANGE1	361
#define	DARKORANGE2	362
#define	DARKORANGE3	363
#define	DARKORANGE4	364
#define	CORAL1	365
#define	CORAL2	366
#define	CORAL3	367
#define	CORAL4	368
#define	TOMATO1	369
#define	TOMATO2	370
#define	TOMATO3	371
#define	TOMATO4	372
#define	ORANGERED1	373
#define	ORANGERED2	374
#define	ORANGERED3	375
#define	ORANGERED4	376
#define	RED1	377
#define	RED2	378
#define	RED3	379
#define	RED4	380
#define	DEEPPINK1	381
#define	DEEPPINK2	382
#define	DEEPPINK3	383
#define	DEEPPINK4	384
#define	HOTPINK1	385
#define	HOTPINK2	386
#define	HOTPINK3	387
#define	HOTPINK4	388
#define	PINK1	389
#define	PINK2	390
#define	PINK3	391
#define	PINK4	392
#define	LIGHTPINK1	393
#define	LIGHTPINK2	394
#define	LIGHTPINK3	395
#define	LIGHTPINK4	396
#define	PALEVIOLETRED1	397
#define	PALEVIOLETRED2	398
#define	PALEVIOLETRED3	399
#define	PALEVIOLETRED4	400
#define	MAROON1	401
#define	MAROON2	402
#define	MAROON3	403
#define	MAROON4	404
#define	VIOLETRED1	405
#define	VIOLETRED2	406
#define	VIOLETRED3	407
#define	VIOLETRED4	408
#define	MAGENTA1	409
#define	MAGENTA2	410
#define	MAGENTA3	411
#define	MAGENTA4	412
#define	ORCHID1	413
#define	ORCHID2	414
#define	ORCHID3	415
#define	ORCHID4	416
#define	PLUM1	417
#define	PLUM2	418
#define	PLUM3	419
#define	PLUM4	420
#define	MEDIUMORCHID1	421
#define	MEDIUMORCHID2	422
#define	MEDIUMORCHID3	423
#define	MEDIUMORCHID4	424
#define	DARKORCHID1	425
#define	DARKORCHID2	426
#define	DARKORCHID3	427
#define	DARKORCHID4	428
#define	PURPLE1	429
#define	PURPLE2	430
#define	PURPLE3	431
#define	PURPLE4	432
#define	MEDIUMPURPLE1	433
#define	MEDIUMPURPLE2	434
#define	MEDIUMPURPLE3	435
#define	MEDIUMPURPLE4	436
#define	THISTLE1	437
#define	THISTLE2	438
#define	THISTLE3	439
#define	THISTLE4	440
#define	GREY0	441
#define	GREY1	442
#define	GREY2	443
#define	GREY3	444
#define	GREY4	445
#define	GREY5	446
#define	GREY6	447
#define	GREY7	448
#define	GREY8	449
#define	GREY9	450
#define	GREY10	451
#define	GREY11	452
#define	GREY12	453
#define	GREY13	454
#define	GREY14	455
#define	GREY15	456
#define	GREY16	457
#define	GREY17	458
#define	GREY18	459
#define	GREY19	460
#define	GREY20	461
#define	GREY21	462
#define	GREY22	463
#define	GREY23	464
#define	GREY24	465
#define	GREY25	466
#define	GREY26	467
#define	GREY27	468
#define	GREY28	469
#define	GREY29	470
#define	GREY30	471
#define	GREY31	472
#define	GREY32	473
#define	GREY33	474
#define	GREY34	475
#define	GREY35	476
#define	GREY36	477
#define	GREY37	478
#define	GREY38	479
#define	GREY39	480
#define	GREY40	481
#define	GREY41	482
#define	GREY42	483
#define	GREY43	484
#define	GREY44	485
#define	GREY45	486
#define	GREY46	487
#define	GREY47	488
#define	GREY48	489
#define	GREY49	490
#define	GREY50	491
#define	GREY51	492
#define	GREY52	493
#define	GREY53	494
#define	GREY54	495
#define	GREY55	496
#define	GREY56	497
#define	GREY57	498
#define	GREY58	499
#define	GREY59	500
#define	GREY60	501
#define	GREY61	502
#define	GREY62	503
#define	GREY63	504
#define	GREY64	505
#define	GREY65	506
#define	GREY66	507
#define	GREY67	508
#define	GREY68	509
#define	GREY69	510
#define	GREY70	511
#define	GREY71	512
#define	GREY72	513
#define	GREY73	514
#define	GREY74	515
#define	GREY75	516
#define	GREY76	517
#define	GREY77	518
#define	GREY78	519
#define	GREY79	520
#define	GREY80	521
#define	GREY81	522
#define	GREY82	523
#define	GREY83	524
#define	GREY84	525
#define	GREY85	526
#define	GREY86	527
#define	GREY87	528
#define	GREY88	529
#define	GREY89	530
#define	GREY90	531
#define	GREY91	532
#define	GREY92	533
#define	GREY93	534
#define	GREY94	535
#define	GREY95	536
#define	GREY96	537
#define	GREY97	538
#define	GREY98	539
#define	GREY99	540
#define	GREY100	541

//const int colorbar[]={RED4,RED,ORANGERED,ORANGE,YELLOW,LAWNGREEN,GREEN,SPRINGGREEN,MEDIUMSPRINGGREEN,TURQUOISE3,DEEPSKYBLUE3,BLUE,NAVYBLUE};
const int colorbar[]={RED,ORANGERED,ORANGE,YELLOW,LAWNGREEN,GREEN,SPRINGGREEN,MEDIUMSPRINGGREEN,TURQUOISE3,DEEPSKYBLUE3,BLUE};
//const int colorbar[]={RED,ORANGERED,ORANGE,YELLOW,RED,YELLOW,GREEN,BLUE,TURQUOISE3,DEEPSKYBLUE3,BLUE};


#endif

