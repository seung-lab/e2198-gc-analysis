function all_types=cell_info_typedef()
load('clusters');

all_types=[
    
%% GC
struct('name','1ni','annotation','','cells',clusters('1ni')'); % 20164
% gallery id: 20164

struct('name','1no','annotation','','cells',clusters('1no')'); % 20092
% gallery id: 20092

% M1
struct('name','1ws','annotation','M1','cells',clusters('1ws')'); % 20203
% gallery id: 20203

% sOffalpha
struct('name','1wt','annotation','sOff\alpha','cells',clusters('1wt')'); % 17109
% gallery id: 17109

% F-mini Off
struct('name','2an','annotation','F-mini^{Off}','cells',clusters('2an')');
% gallery id: 15066

% F-midi Off, J
struct('name','2aw','annotation','F-midi^{Off},J','cells',clusters('2aw')'); 
% gallery id: 17028

struct('name','2i','annotation','','cells',clusters('2i')');
% gallery id: 20234

struct('name','2o','annotation','','cells',clusters('2o')');
% gallery id: 10005

struct('name','25','annotation','','cells',clusters('25')');
% gallery id: 17132

struct('name','27','annotation','','cells',clusters('27')'); 
% gallery id: 26065

struct('name','28','annotation','','cells',clusters('28')'); 
% gallery id: 20155

struct('name','3i','annotation','','cells',clusters('3i')');
% gallery id: 17135

struct('name','3o','annotation','','cells',clusters('3o')');
% gallery id: 17037

% On-Off DS caudal
struct('name','37c','annotation','On-Off DS','cells',clusters('37c')'); 
% gallery id: 20210

% On-Off DS dorsal
struct('name','37d','annotation','On-Off DS','cells',clusters('37d')'); 
% gallery id: 20125

% On-Off DS rostral
struct('name','37r','annotation','On-Off DS','cells',clusters('37r')'); 
% gallery id: 20213

% On-Off DS ventral
struct('name','37v','annotation','On-Off DS','cells',clusters('37v')'); 
% gallery id: 20233

struct('name','4i','annotation','','cells',clusters('4i')'); 
% gallery id: 20170

% mini-tOFFalpha
struct('name','4on','annotation','','cells',clusters('4on')'); 
% gallery id: 20230

% tOffalpha
struct('name','4ow','annotation','tOff\alpha','cells',clusters('4ow')'); 
% gallery id: 20156
 
struct('name','5si','annotation','HD1','cells',clusters('5si')'); 
%  gallery id: 20183

struct('name','5so','annotation','HD2','cells',clusters('5so')');
% gallery id: 20223

struct('name','5ti','annotation','UHD','cells',clusters('5ti')');
% gallery id: 20184

struct('name','5to','annotation','','cells',clusters('5to')'); 
% gallery id: 20240

% W3 
struct('name','51','annotation','W3/LED','cells',clusters('51')');
% gallery id: 20212

% F-mini On
struct('name','63','annotation','F-mini^{On}','cells',clusters('63')');
% gallery id: 20178

% mini-tONalpha
struct('name','6sn','annotation','','cells',clusters('6sn')'); 
% gallery id: 20198

% tOnalpha
struct('name','6sw','annotation','tOn\alpha','cells',clusters('6sw')'); 
% gallery id: 20222

% F-midi On 
struct('name','6t','annotation','F-midi^{On}','cells',clusters('6t')'); 
% gallery id: 20232

% On DS
struct('name','7id','annotation','sOn DS','cells',clusters('7id')'); 
% gallery id: 26078

% On DS
struct('name','7ir','annotation','sOn DS','cells',clusters('7ir')'); 
% gallery id: 26002

% On DS
struct('name','7iv','annotation','sOn DS','cells',clusters('7iv')'); 
% gallery id: 26077 

% On DS
struct('name','7o','annotation','tOn DS','cells',clusters('7o')'); 
% gallery id: 20239

struct('name','72','annotation','','cells',clusters('72')'); 
% gallery id: 20221

struct('name','73','annotation','OND','cells',clusters('73')'); 
% gallery id: 20100

struct('name','8n','annotation','','cells',clusters('8n')'); 
% gallery id: 20126

% sOnalpha/M4
struct('name','8w','annotation','sON\alpha/M4','cells',clusters('8w')'); 
% gallery id: 17111

struct('name','81i','annotation','','cells',clusters('81i')'); 
% gallery id: 20158

struct('name','81o','annotation','','cells',clusters('81o')'); 
% gallery id: 20069

struct('name','82n','annotation','','cells',clusters('82n')'); 
% gallery id: 20161

struct('name','82wi','annotation','','cells',clusters('82wi')'); 
% gallery id: 30001

% regular bilayer 
struct('name','82wo','annotation','','cells',clusters('82wo')');
% gallery id: 20118

struct('name','85','annotation','','cells',clusters('85')'); 
% gallery id: 20200

% M5
struct('name','9n','annotation','','cells',clusters('9n')'); 
% gallery id: 20112

% M2
struct('name','9w','annotation','M2','cells',clusters('9w')'); 
% gallery id: 20228

% M3
struct('name','91','annotation','','cells',clusters('91')'); 
% gallery id: 20218

struct('name','915','annotation','','cells',clusters('915')'); 
% gallery id: 17009

struct('name','weirdos','annotation','','cells',[17134 20248 26069 26093]);

struct('name','cutoffs','annotation','','cells',[10015 17051 17238 26081 26153 26030 26105 26114 26176 26107 26108]);
    
%% BC
% 1 
struct('name','bc1','annotation','NK3R+ & Syt2-','cells',[60008 60019 60026 60027 60032 60052 60055 60078 ...
                     60079 60099 60105 60109 60110 60111 60114 60118 60129 60132 60139 60142 ...
                     60147 60150 60158 60161 60162 60164 60171 60177 60184 60187 60188 60189 ...
                     60194 60195 60196 60203 60204 60212 60213 60216 60218]);
% 2 
struct('name','bc2','annotation','NK3R+ & Syt2+,recoverin,Neto1,Cdh8','cells',[60080 60001 60002 60003 60004 60006 60009 60010 ...
                     60011 60012 60013 60021 60022 60023 60025 60029 60037 60038 60039 60040 ...
                     60041 60042 60043 60044 60046 60097 60101 60102 60103 60104 60112 60120 ...
                     60124 60130 60133 60135 60138 60140 60141 60149 60157 60159 60167 60169 ...
                     60182 60208 60209 60210 60214 60217 60219 60221 60223 60224 60226]);
% 3 
struct('name','bc3a','annotation','HCN4','cells',[60028 60030 60048 60049 60050 60056 60059 60066 ...
                     60068 60072 60075 60076 60085 60088 60092 60093 60108 60115 60123 60134 ...
                     60136 60137 60145 60146 60166 60172 60176 60181]);
% 4 
struct('name','bc3b','annotation','PKARIIbeta','cells',[60024 60045 60053 60057 60063 60069 60073 60074 ...
                     60077 60082 60084 60087 60090 60094 60096 60098 60106 60107 60122 60131 ...
                     60148 60151 60153 60154 60156 60165 60173 60178 60179 60190 60201 60211 ...
                     60215 60228]);
% 5 
struct('name','bc4','annotation','Csen','cells',[60047 60058 60067 60083 60086 60089 60091 60095 ...
                     60100 60116 60117 60119 60121 60125 60127 60143 60144 60160 60163 60168 ...
                     60171 60174 60175 60183 60185 60186 60191 60197 60199 60202 60206 60220 ...
                     60222]);
% 6 
struct('name','bc5t','annotation','','cells',[60543 60452 60366 60371 60387 60395 60399 60403 ...
                      60408 60411 60436 60445 60448 60465 60469 60481 60502 60503 60507 60527 ...
                      60533 60535 60538 60054 60155]);
% 7 
struct('name','bc5o','annotation','5-HT3R-EGFP,Kcng4,Cdh9','cells',[60364 60390 60405 60406 60420 60432 60434 60440 ...
                      60447 60453 60459 60466 60473 60484 60490 60495 60513 60534 60540 60556 ...
                      60559 60608 60612 60018 60061 60064 60071 60532]);
% 8 
struct('name','bc5i','annotation','5-HT3R-EGFP,Kcng4,Cdh9','cells',[60360 60033 60374 60380 60383 60386 60388 60389 ...
                     60404 60410 60414 60415 60421 60439 60442 60450 60458 60460 60462 60478 ...
                     60488 60491 60497 60498 60504 60505 60510 60514 60519 60522 60523 60528 ...
                     60541 60542 60020 60615 60617 60619 60620 60621 60618]);
% 9 
struct('name','xbc','annotation','','cells',[60355 60449 60379 60413 60430 60455 60493 60501 ...
                     60517 60539 60547 60550]);
% 10 
struct('name','bc6','annotation','Syt2','cells',[60489 60422 60356 60060 60423 60394 60425 60398 ...
                     60401 60426 60407 60428 60443 60444 60476 60477 60483 60486 60496 60499 ...
                     60512 60516 60526 60548 60549 60552 60554 60561 60611 60614 60467 60358 ...
                     60370 60375 60381 60416 60419 60431 60441 60451 60464 60487 60530 60537 ...
                     60036 60363 60409]);
% 11 
struct('name','bc7','annotation','Gus-GFP','cells',[60361 60373 60376 60377 60051 60393 60396 60397 ...
                     60412 60418 60429 60435 60437 60446 60454 60470 60480 60485 60492 60508 ...
                     60509 60525 60529 60531 60546 60553 60562 60354 60016]);
% 12 
struct('name','bc8/9','annotation','/Clm1','cells',[60433 60482 60500 60578 60368 60402 60417 60438 ...
                       60457 60461]);
% 13 
struct('name','rbc','annotation','','cells',[60463 60536 60357 60359 60365 60369 60372 60378 ...
                     60382 60384 60392 60400 60456 60468 60471 60472 60474 60475 60479 60506 ...
                     60515 60518 60520 60521 60544 60551 60555 60557 60558 60560 60563 60564 ...
                     60565 60566 60567 60568 60569 60570 60571 60572 60573 60574 60575 60576 ...
                     60577 60580 60581 60582 60583 60584 60585 60586 60587 60588 60589 60590 ...
                     60591 60592 60593 60594 60595 60596 60597 60598 60599 60600 60601 60602 ...
                     60603 60604 60605 60606 60607 60609 60610 60613 60017 60031 60427]);
                     
%% AC
% 1 On SAC
struct('name','ON SAC','annotation','','cells', ...
           [70244 70243 70242 70241 70240 70239 70238 70237 70236 70235 70234 70233 ...
            70232 70231 70230 70229 70228 70227 70225 70224 70223 70222 70221 70220 ...
            70219 70218 70217 70216 70215 70214 70213 70212 70211 70209 70208 70207 ...
            70206 70205 70204 70203 70202 70201 70200 70199 70198 70197 70196 70195 ...
            70194 70193 70192 70191 70189 70188 70187 70186 70185 70184 70183 70182 ...
            70181 70180 70179 70178 70176 70174 70172 70171 70170 70169 70168 70164 ...
            70163 70162 70161 70029 70028 20204 20196 20062 20044 20032 20030 20025 ...
            26007 26009 26010 26011 26012 26013 26014 26015 26016 26017 26139 26143 ...
            26169 26174 26180 26183 26184 26186 26197]);

% 2 Off SAC
struct('name','OFF SAC','annotation','','cells', ...
           [70167 70158 70156 70155 70154 70152 70151 70150 70149 70148 70147 70146 ...
            70145 70141 70138 70137 70134 70133 70131 70130 70129 70128 70127 70126 ...
            70125 70124 70123 70122 70121 70120 70119 70118 70117 70116 70115 70114 ...
            70113 70112 70111 70110 70109 70108 70106 70102 70100 70099 70096 70095 ...
            70093 70090 70089 70088 70087 70086 70085 70084 70083 70082 70081 70080 ...
            70079 70077 70076 70068 70066 70050 70048 70035 70034 70033 70032 70031 ...
            70030 70027 70026 70025 70024 70023 70016 70014]);

% bifur
% 40001 40002 41004 41005 41053 41041 41045 41065 41082 41098 41096

];

end