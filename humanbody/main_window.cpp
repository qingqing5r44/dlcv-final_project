#include<iostream>
#include <stdio.h>   //popen()
//#include <string.h>  //memset()
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
//#include "SparsePCAForShapeGrad.h"
#include <omp.h>
#include "faceBodyModel.h"
#include "SparseLocalizedForConnectMap.h"
#include "edit.h"
#include <Python.h>
#include "subdivision.h"
#include <math.h>
#include<time.h>
#include <QFileDialog>
#include <QtCore/QTextCodec>
#include <QtCore/qdebug.h>
#include <windows.h>
#include <GL/glut.h>

#include "main_window.h"
#pragma execution_character_set("utf-8")
#pragma comment(lib,"opengl32.lib")
#pragma comment(lib,"glu32.lib")
#pragma comment(lib,"glut32.lib")

using namespace std;
typedef  OpenMesh::PolyMesh_ArrayKernelT<> PGMesh;
void preSubdivisionDatasets();
void processingDETrainingSets1();
struct Mesh
{
	Eigen::MatrixXd V, U, C, VC;
	Eigen::MatrixXi F;
}low, high, scene;
std::vector<integer> selected_vid_vec;
decimal color_fac = 40;
integer switch_ld = 0;
int faceSet[1719] = {
	3902,
	3903,
	3906,
	3907,
	3918,
	3919,
	3930,
	3931,
	3936,
	3937,
	3940,
	3941,
	3942,
	3943,
	3960,
	3961,
	3966,
	3967,
	3986,
	3987,
	4110,
	4111,
	4112,
	4113,
	4114,
	4115,
	4116,
	4117,
	4118,
	4119,
	4120,
	4121,
	4122,
	4123,
	4124,
	4125,
	4126,
	4127,
	4128,
	4129,
	4130,
	4131,
	4132,
	4133,
	4134,
	4135,
	4136,
	4137,
	4138,
	4139,
	4140,
	4141,
	4142,
	4143,
	4144,
	4145,
	4146,
	4147,
	4148,
	4149,
	4286,
	4287,
	4292,
	4293,
	4294,
	4295,
	4296,
	4297,
	4298,
	4299,
	4328,
	4329,
	4330,
	4331,
	4332,
	4333,
	4334,
	4335,
	4336,
	4337,
	4338,
	4339,
	4342,
	4343,
	4344,
	4345,
	4408,
	4409,
	4548,
	4549,
	4550,
	4551,
	4552,
	4553,
	4554,
	4555,
	4556,
	4557,
	4558,
	4559,
	4560,
	4561,
	4562,
	4563,
	4564,
	4565,
	4566,
	4567,
	4568,
	4569,
	4570,
	4571,
	4572,
	4573,
	4574,
	4575,
	4576,
	4577,
	4578,
	4579,
	4580,
	4581,
	4582,
	4583,
	4584,
	4585,
	4586,
	4587,
	4588,
	4589,
	4590,
	4591,
	4592,
	4593,
	4594,
	4595,
	4596,
	4597,
	4598,
	4599,
	4600,
	4601,
	4602,
	4603,
	4700,
	4701,
	4704,
	4705,
	4716,
	4717,
	4728,
	4729,
	4734,
	4735,
	4738,
	4739,
	4740,
	4741,
	4758,
	4759,
	4764,
	4765,
	4784,
	4785,
	4908,
	4909,
	4910,
	4911,
	4912,
	4913,
	4914,
	4915,
	4916,
	4917,
	4918,
	4919,
	4920,
	4921,
	4922,
	4923,
	4924,
	4925,
	4926,
	4927,
	4928,
	4929,
	4930,
	4931,
	4932,
	4933,
	4934,
	4935,
	4936,
	4937,
	4938,
	4939,
	4940,
	4941,
	4942,
	4943,
	4944,
	4945,
	4946,
	4947,
	5084,
	5085,
	5090,
	5091,
	5092,
	5093,
	5094,
	5095,
	5096,
	5097,
	5126,
	5127,
	5128,
	5129,
	5130,
	5131,
	5132,
	5133,
	5134,
	5135,
	5136,
	5137,
	5140,
	5141,
	5142,
	5143,
	5206,
	5207,
	5320,
	5321,
	5322,
	5323,
	5324,
	5325,
	5326,
	5327,
	5328,
	5329,
	5330,
	5331,
	5332,
	5333,
	5334,
	5335,
	5336,
	5337,
	5338,
	5339,
	5340,
	5341,
	5342,
	5343,
	5344,
	5345,
	5346,
	5347,
	5348,
	5349,
	5350,
	5351,
	5352,
	5353,
	5354,
	5355,
	5356,
	5357,
	5358,
	5359,
	5360,
	5361,
	5362,
	5363,
	5364,
	5365,
	5366,
	5367,
	5368,
	5369,
	5370,
	5371,
	5372,
	5373,
	5374,
	5375,
	5646,
	5647,
	5648,
	5649,
	5650,
	5651,
	5652,
	5653,
	5654,
	5655,
	5658,
	5659,
	5660,
	5661,
	5663,
	5664,
	5665,
	5666,
	5667,
	5668,
	5669,
	5670,
	5671,
	5672,
	5673,
	5674,
	5675,
	5676,
	5677,
	5678,
	5679,
	5680,
	5681,
	5682,
	5683,
	5684,
	5685,
	5686,
	5687,
	5688,
	5689,
	5690,
	5691,
	5692,
	5693,
	5696,
	5697,
	5698,
	5699,
	5701,
	5702,
	5703,
	5704,
	5705,
	5706,
	5707,
	5708,
	5709,
	5710,
	5711,
	5712,
	5713,
	5714,
	5715,
	5716,
	5717,
	5718,
	5719,
	5720,
	5721,
	5722,
	5723,
	5724,
	5725,
	5726,
	5727,
	5728,
	5729,
	5732,
	5733,
	5734,
	5735,
	5736,
	5737,
	5738,
	5739,
	5740,
	5741,
	5742,
	5743,
	5744,
	5745,
	5746,
	5747,
	5748,
	5749,
	5750,
	5751,
	5752,
	5753,
	5754,
	5755,
	5756,
	5757,
	5758,
	5759,
	5760,
	5761,
	5762,
	5763,
	5764,
	5765,
	5766,
	5767,
	5768,
	5769,
	5770,
	5771,
	5772,
	5773,
	5774,
	5775,
	5776,
	5777,
	5778,
	5779,
	5780,
	5781,
	5782,
	5783,
	5784,
	5785,
	5786,
	5787,
	5788,
	5789,
	5790,
	5791,
	5792,
	5793,
	5794,
	5795,
	5796,
	5797,
	5798,
	5799,
	5800,
	5801,
	5802,
	5803,
	5804,
	5805,
	5806,
	5807,
	5808,
	5809,
	5810,
	5811,
	5812,
	5813,
	5814,
	5815,
	5816,
	5817,
	5818,
	5819,
	5820,
	5821,
	5822,
	5823,
	5824,
	5825,
	5826,
	5827,
	5828,
	5829,
	5830,
	5831,
	5832,
	5833,
	5834,
	5835,
	5836,
	5837,
	5838,
	5839,
	5840,
	5841,
	5842,
	5843,
	5844,
	5845,
	5846,
	5847,
	5848,
	5849,
	5850,
	5851,
	5852,
	5853,
	5854,
	5855,
	5856,
	5857,
	5858,
	5859,
	5860,
	5861,
	5862,
	5863,
	5864,
	5865,
	5866,
	5867,
	5868,
	5869,
	5870,
	5871,
	5872,
	5873,
	5874,
	5875,
	5876,
	5877,
	5878,
	5879,
	5880,
	5881,
	5882,
	5883,
	5884,
	5885,
	5886,
	5887,
	5888,
	5889,
	5890,
	5891,
	5892,
	5893,
	5894,
	5895,
	5896,
	5897,
	5898,
	5899,
	5900,
	5901,
	5902,
	5903,
	5904,
	5905,
	5906,
	5907,
	5908,
	5909,
	5910,
	5911,
	5912,
	5913,
	5914,
	5915,
	5916,
	5917,
	5918,
	5919,
	5920,
	5921,
	5922,
	5923,
	5924,
	5925,
	5926,
	5927,
	5928,
	5929,
	5930,
	5931,
	5932,
	5933,
	5934,
	5935,
	5936,
	5937,
	5938,
	5939,
	5940,
	5941,
	5942,
	5943,
	5944,
	5945,
	5946,
	5947,
	5948,
	5949,
	5950,
	5951,
	5952,
	5953,
	5954,
	5955,
	5956,
	5957,
	5958,
	5959,
	5960,
	5961,
	5962,
	5963,
	5964,
	5965,
	5966,
	5967,
	5968,
	5969,
	5970,
	5971,
	5972,
	5973,
	5974,
	5975,
	5976,
	5977,
	5978,
	5979,
	5980,
	5981,
	5982,
	5983,
	5984,
	5985,
	5986,
	5987,
	5988,
	5989,
	5990,
	5991,
	5992,
	5993,
	5994,
	5995,
	5996,
	5997,
	5998,
	5999,
	6000,
	6001,
	6002,
	6003,
	6004,
	6005,
	6006,
	6007,
	6008,
	6009,
	6010,
	6011,
	6012,
	6013,
	6014,
	6015,
	6016,
	6017,
	6018,
	6019,
	6020,
	6021,
	6022,
	6023,
	6024,
	6025,
	6026,
	6027,
	6028,
	6029,
	6030,
	6031,
	6032,
	6033,
	6034,
	6035,
	6036,
	6037,
	6038,
	6039,
	6040,
	6041,
	6042,
	6043,
	6044,
	6045,
	6046,
	6047,
	6048,
	6049,
	6050,
	6051,
	6052,
	6053,
	6054,
	6055,
	6056,
	6057,
	6058,
	6059,
	6060,
	6061,
	6062,
	6063,
	6064,
	6065,
	6066,
	6067,
	6068,
	6069,
	6070,
	6071,
	168,
	169,
	170,
	171,
	172,
	173,
	174,
	175,
	176,
	177,
	178,
	179,
	180,
	181,
	182,
	183,
	184,
	185,
	186,
	187,
	188,
	189,
	190,
	191,
	192,
	193,
	194,
	195,
	196,
	197,
	198,
	199,
	200,
	201,
	202,
	203,
	204,
	205,
	206,
	207,
	208,
	209,
	210,
	211,
	212,
	213,
	214,
	215,
	216,
	217,
	218,
	219,
	220,
	221,
	222,
	223,
	224,
	225,
	226,
	227,
	228,
	229,
	230,
	231,
	232,
	233,
	234,
	235,
	236,
	237,
	238,
	239,
	240,
	241,
	242,
	243,
	244,
	245,
	246,
	247,
	248,
	249,
	250,
	251,
	252,
	253,
	254,
	255,
	256,
	257,
	258,
	259,
	260,
	261,
	262,
	263,
	264,
	265,
	266,
	267,
	268,
	269,
	270,
	271,
	272,
	273,
	274,
	275,
	276,
	277,
	278,
	279,
	280,
	281,
	282,
	283,
	284,
	285,
	286,
	287,
	288,
	289,
	290,
	291,
	292,
	293,
	294,
	295,
	296,
	297,
	298,
	299,
	300,
	301,
	302,
	303,
	304,
	305,
	306,
	307,
	308,
	309,
	310,
	311,
	312,
	313,
	314,
	315,
	316,
	317,
	318,
	319,
	320,
	321,
	322,
	323,
	324,
	325,
	326,
	327,
	328,
	329,
	330,
	331,
	332,
	333,
	334,
	335,
	336,
	337,
	338,
	339,
	340,
	341,
	342,
	343,
	344,
	345,
	346,
	347,
	348,
	349,
	350,
	351,
	352,
	353,
	354,
	355,
	356,
	357,
	358,
	359,
	360,
	361,
	362,
	363,
	364,
	365,
	366,
	367,
	368,
	369,
	370,
	371,
	372,
	373,
	374,
	375,
	376,
	377,
	378,
	379,
	380,
	381,
	382,
	383,
	384,
	385,
	386,
	387,
	388,
	389,
	390,
	391,
	392,
	393,
	394,
	395,
	396,
	397,
	398,
	399,
	400,
	401,
	402,
	403,
	404,
	405,
	406,
	407,
	408,
	409,
	410,
	411,
	412,
	413,
	414,
	415,
	416,
	417,
	418,
	419,
	420,
	421,
	422,
	423,
	424,
	425,
	426,
	427,
	428,
	429,
	430,
	431,
	432,
	433,
	434,
	435,
	436,
	437,
	438,
	439,
	440,
	441,
	442,
	443,
	444,
	445,
	446,
	447,
	448,
	449,
	450,
	451,
	452,
	453,
	454,
	455,
	456,
	457,
	458,
	459,
	460,
	461,
	462,
	463,
	464,
	465,
	466,
	467,
	468,
	469,
	470,
	471,
	472,
	473,
	474,
	475,
	476,
	477,
	478,
	479,
	480,
	481,
	482,
	483,
	484,
	485,
	486,
	487,
	488,
	489,
	490,
	491,
	492,
	493,
	494,
	495,
	496,
	497,
	498,
	499,
	500,
	501,
	502,
	503,
	504,
	505,
	506,
	507,
	508,
	509,
	510,
	511,
	512,
	513,
	514,
	515,
	516,
	517,
	518,
	519,
	520,
	521,
	522,
	523,
	524,
	525,
	526,
	527,
	528,
	529,
	530,
	531,
	532,
	533,
	534,
	535,
	536,
	537,
	538,
	539,
	540,
	541,
	542,
	543,
	544,
	545,
	546,
	547,
	548,
	549,
	550,
	551,
	552,
	553,
	554,
	555,
	556,
	557,
	558,
	559,
	560,
	561,
	562,
	563,
	564,
	565,
	566,
	567,
	568,
	569,
	570,
	571,
	572,
	573,
	574,
	575,
	576,
	577,
	578,
	579,
	580,
	581,
	582,
	583,
	19072,
	19073,
	19074,
	19075,
	19076,
	19077,
	19078,
	19079,
	19080,
	19081,
	19082,
	19083,
	19084,
	19085,
	19086,
	19087,
	19088,
	19089,
	19090,
	19091,
	19092,
	19093,
	19094,
	19095,
	19096,
	19097,
	19098,
	19099,
	19100,
	19101,
	19102,
	19103,
	19104,
	19105,
	19106,
	19107,
	19108,
	19109,
	19110,
	19111,
	19112,
	19113,
	19114,
	19115,
	19116,
	19117,
	19118,
	19119,
	19120,
	19121,
	19122,
	19123,
	19124,
	19125,
	19126,
	19127,
	19128,
	19129,
	19130,
	19131,
	19132,
	19133,
	19134,
	19135,
	19136,
	19137,
	19138,
	19139,
	19140,
	19141,
	19142,
	19143,
	19144,
	19145,
	19146,
	19147,
	19148,
	19149,
	19150,
	19151,
	19152,
	19153,
	19154,
	19155,
	19156,
	19157,
	19158,
	19159,
	19160,
	19161,
	19162,
	19163,
	19164,
	19165,
	19166,
	19167,
	19168,
	19169,
	19170,
	19171,
	19172,
	19173,
	19174,
	19175,
	19176,
	19177,
	19178,
	19179,
	19180,
	19181,
	19182,
	19183,
	19184,
	19185,
	19186,
	19187,
	19188,
	19189,
	19190,
	19191,
	19192,
	19193,
	19194,
	19195,
	19196,
	19197,
	19198,
	19199,
	19250,
	19251,
	19255,
	19256,
	19257,
	19258,
	19259,
	19260,
	19261,
	19262,
	19263,
	19264,
	19265,
	19266,
	19267,
	19268,
	19269,
	19270,
	19271,
	19325,
	19328,
	19329,
	19330,
	19331,
	19332,
	19333,
	19334,
	19335,
	19340,
	19341,
	19342,
	19343,
	19348,
	19349,
	19354,
	19355,
	19356,
	19357,
	19358,
	19359,
	19364,
	19372,
	19373,
	19374,
	19375,
	19376,
	19392,
	19393,
	19394,
	19395,
	19396,
	19397,
	19398,
	19399,
	19402,
	19403,
	19446,
	19447,
	19448,
	19449,
	19450,
	19451,
	19483,
	19487,
	19501,
	19502,
	19503,
	19504,
	19505,
	19506,
	19507,
	19508,
	19509,
	19510,
	19511,
	19512,
	19513,
	19514,
	19515,
	19530,
	19531,
	19532,
	19534,
	19535,
	19536,
	19537,
	19538,
	19539,
	19540,
	19541,
	19542,
	19543,
	19560,
	19562,
	19563,
	20514,
	20515,
	20522,
	20523,
	20529,
	20530,
	20531,
	20546,
	20547,
	20548,
	20549,
	20550,
	20551,
	20552,
	20553,
	20554,
	20555,
	20556,
	20557,
	20558,
	20559,
	20560,
	20561,
	21122,
	21123,
	21130,
	21131,
	21137,
	21138,
	21139,
	21154,
	21155,
	21156,
	21157,
	21158,
	21159,
	21160,
	21161,
	21162,
	21163,
	21164,
	21165,
	21166,
	21167,
	21168,
	21169,
	22104,
	22105,
	22106,
	22107,
	22108,
	22109,
	22113,
	22114,
	22115,
	22120,
	22121,
	22128,
	22129,
	22166,
	22167,
	22168,
	22169,
	22170,
	22171,
	22172,
	22173,
	22174,
	22175,
	22176,
	22177,
	22178,
	22179,
	22180,
	22181,
	22202,
	22203,
	22234,
	22235,
	22236,
	22237,
	22238,
	22239,
	22240,
	22241,
	22242,
	22243,
	22244,
	22245,
	22246,
	22247,
	22248,
	22249,
	22250,
	22251,
	22252,
	22253,
	22254,
	22255,
	22256,
	22257,
	22258,
	22259,
	22263,
	22264,
	22265,
	22266,
	22267,
	22268,
	22269,
	22270,
	22271,
	22272,
	22273,
	22274,
	22275,
	23040,
	23041,
	23044,
	23045,
	23050,
	23051,
	23057,
	23064,
	23065,
	23075,
	23102,
	23103,
	23104,
	23105,
	23106,
	23107,
	23108,
	23109,
	23110,
	23111,
	23112,
	23113,
	23114,
	23115,
	23116,
	23117,
	23138,
	23139,
	23170,
	23171,
	23172,
	23173,
	23174,
	23175,
	23176,
	23177,
	23178,
	23179,
	23180,
	23181,
	23182,
	23183,
	23184,
	23185,
	23186,
	23187,
	23188,
	23189,
	23190,
	23191,
	23192,
	23193,
	23194,
	23195,
	23196,
	23197,
	23198,
	23199,
	23200,
	23201,
	23202,
	23203,
	23204,
	23205,
	23206,
	23207,
	23208,
	23209,
	23210,
	23211,
	23997,
	24132,
	24133,
	24134,
	24135,
	24136,
	24137,
	24138,
	24139,
	24140,
	24141,
	24142,
	24143,
	24144,
	24145,
	24146,
	24147,
	24168,
	24169,
	24172,
	24173,
	24178,
	24179,
	24180,
	24181,
	24182,
	24183,
	24232,
	24274,
	24275,
	24278,
	24279,
	24280,
	24281,
	24282,
	24283,
	24284,
	24285,
	24286,
	24287,
	24288,
	24289,
	24290,
	24291,
	24292,
	24293,
	24294,
	24295,
	24296,
	24297,
	24298,
	24299,
	24300,
	24301,
	24302,
	24303,
	24304,
	24305,
	24306,
	24307,
	24344,
	24345,
	24348,
	24349,
	24350,
	24351,
	24396,
	24484,
	24485,
	24486,
	24487,
	24488,
	24489,
	24490,
	24491,
	24492,
	24493,
	24494,
	24495,
	24496,
	24497,
	24498,
	24499,
	24508,
	24509,
	24510,
	24511,
	24518,
	24520,
	24521,
	24524,
	24525,
	24530,
	24531,
	24532,
	24533,
	24534,
	24535,
	24584,
	24585,
	24597,
	24612,
	24613,
	24614,
	24615,
	24626,
	24627,
	24630,
	24631,
	24632,
	24633,
	24634,
	24635,
	24636,
	24637,
	24638,
	24639,
	24640,
	24641,
	24642,
	24643,
	24644,
	24645,
	24646,
	24647,
	24648,
	24649,
	24650,
	24651,
	24652,
	24653,
	24654,
	24655,
	24656,
	24657,
	24658,
	24659,
	24696,
	24697,
	24712,
	24713,
	24736,
	24737,
	24738,
	24739,
	24748,
	24749,
	24752,
	24753,
	24766,
	24767,
	24803,
	24809,
	24810,
	24811,
	24812,
	24813,
	24815,
	24816,
	24817,
	24822,
	24823,
	24836,
	24837,
	24892,
	24893,
	24912,
	24913,
	24918,
	24919,
	24921,
	24922,
	24923,
	24934,
	24935,
	24983,
	24984,
	24985,
	24991,
	24993,
	24994,
	24995,
	24998,
	24999,
	25004,
	25005,
	25018,
	25019
};
typedef OpenMesh::Vec3f float3;
std::vector<float3> neutralExpr(11510);					// B_0
std::vector<std::vector<float3>> exprList(46);		// B_i ( 1 <= i <= 46 )
				
													
													
													// fileName: e.g. "shape.bs"


//
Eigen::VectorXd getRLA(Eigen::VectorXd x0, Eigen::VectorXd x1)
{
	int nE = x0.size() / 2;
	Eigen::VectorXd x(2 * nE);
	
	for (int i = 0; i < nE; i++)
	{
		x(2 * i + 0) = x1(2 * i + 0) - x0(2 * i + 0);
		x(2 * i + 1) = x1(2 * i + 1) / x0(2 * i + 1) - 1;
	}
	return x;
}

Eigen::VectorXd getLAfromRLA(Eigen::VectorXd x0, Eigen::VectorXd x)
{
	int nE = x0.size() / 2;
	Eigen::VectorXd x1(nE* 2);

	for (int i = 0; i < nE; i++)
	{
		x1(2 * i + 0) = x(2 * i + 0) + x0(2 * i + 0);
		x1(2 * i + 1) =( x(2 * i + 1)+ 1) * x0(2 * i + 1) ;
	}
	return x1;
}


void loadBlendshape(const std::string& fileName)
{
	FILE* fp;
	fopen_s(&fp, fileName.c_str(), "rb");

	int nShapes = 0, nVerts = 0, nFaces = 0;
	fread(&nShapes, sizeof(int), 1, fp);			// nShape = 46
	fread(&nVerts, sizeof(int), 1, fp);			// nVerts = 11510
	fread(&nFaces, sizeof(int), 1, fp);			// nFaces = 11540

												// Load neutral expression B_0
	fread(&neutralExpr[0], sizeof(float3), nVerts, fp);

	// Load other expressions B_i ( 1 <= i <= 46 )
	for (int exprId = 0; exprId<nShapes; exprId++) {
		std::vector<float3>& expr = exprList[exprId];
		expr.resize(nVerts);
		fread(&expr[0], sizeof(float3), nVerts, fp);
	}

	fclose(fp);
}
void reconstructionFromNComp();
void saveExpressionTensor()
{
	char *meshpath = new char[100];
	//sprintf_s(meshpath, 100, "F:\\scan\\male_without_mouth\\0000.obj");
	sprintf_s(meshpath, 100, "F:\\facewarehouse\\shapebasis.txt");
	ofstream file(meshpath);
	for (int i_th = 0; i_th < 150; i_th++)
	{
		sprintf_s(meshpath, 100, "F:\\facewarehouse\\Tester_%d\\Blendshape\\shape.bs", i_th + 1);
		loadBlendshape(meshpath);
		for (int j = 0; j < neutralExpr.size(); j++)
		{
			file << neutralExpr[j][0] << " " << neutralExpr[j][1]  << " " << neutralExpr[j][2] << " ";
		}
		file << endl;
		
		sprintf_s(meshpath, 100, "F:\\facewarehouse\\%d.txt", i_th);
		ofstream file1(meshpath);
		for (int i = 0; i < exprList.size(); i++)
		{
			for (int j = 0; j < neutralExpr.size(); j++)
			{
				file1 << exprList[i][j][0] - neutralExpr[j][0] << " " << exprList[i][j][1] - neutralExpr[j][1] << " " << exprList[i][j][2] - neutralExpr[j][2] << " ";
			}
			file1 << endl;
		}
		file1.close();
	}
}

void readBodyFeaturePos(string  p_path, vector<Eigen::Vector3d> &body_feature_pos)
{
	ifstream file(p_path);
	for (int i = 0; i < 19; i++)
	{
		Eigen::Vector3d content;
		file >> content(0);
		file >> content(1);
		file >> content(2);
		body_feature_pos[i] = content;
	}
	return;
}

void readBodyFeaturePos1(string  p_path, vector<Eigen::Vector3d> &body_feature_pos)
{
	ifstream file(p_path);
	for (int i = 0; i < 33; i++)
	{
		Eigen::Vector3d content;
		file >> content(0);
		file >> content(1);
		file >> content(2);
		body_feature_pos[i] = content;
	}
	return;
}

void readHandFeaturePos(string  p_path, vector<Eigen::Vector3d> &hand_feature_pos, vector<vector<int>> & hand_feature_idx, vector<vector<double>> & hand_feature_weight)
{
	ifstream file(p_path);
	int count;
	file >> count;
	hand_feature_pos.resize(0);
	hand_feature_idx.resize(0);
	hand_feature_weight.resize(0);
	while (count > 0)
	{
		
		/*hand_feature_idx[i].resize(count);
		hand_feature_weight[i].resize(count);*/
		vector<int> idx(count);
		vector<double> weight(count);
		for (int j = 0; j < count; j++)
		{
			file >> idx[j];
			weight[j] = 1.0 / count;
		}
		Eigen::Vector3d content;
		file >> content(0);
		file >> content(1);
		file >> content(2);
		hand_feature_pos.push_back( content);
		hand_feature_idx.push_back(idx);
		hand_feature_weight.push_back(weight);
		file >> count;
	}
	return;
}

void readFacecontourFeaturePos(string  p_path, vector<Eigen::Vector3d> &hand_feature_pos, vector<vector<int>> & hand_feature_idx, vector<vector<double>> & hand_feature_weight)
{
	ifstream file(p_path);
	int count;
	file >> count;
	while (count != -1)
	{
		vector<int> idx(count);
		vector<double> weight(count);
		for (int j = 0; j < count; j++)
		{
			file >> idx[j];
			weight[j] = 1.0 / count;
		}
		Eigen::Vector3d content;
		file >> content(0);
		file >> content(1);
		file >> content(2);
		hand_feature_pos.push_back(content);
		hand_feature_idx.push_back(idx);
		hand_feature_weight.push_back(weight);
		file >> count;
	}
	return;
}


void readCorespondence(string  p_path, vector<Eigen::Vector3d> &hand_feature_pos, vector<vector<int>> & hand_feature_idx, vector<vector<double>> & hand_feature_weight)
{
	hand_feature_pos.resize(0);
	hand_feature_idx.resize(0);
	hand_feature_weight.resize(0);
	ifstream file(p_path);
	int x = 0;
	int count;
	file >> count;
	while(count >= 0)
	{
		vector<int> idx = { count };
		vector<double> weight = { 1.0 };
		Eigen::Vector3d pos;
		file >> pos[0];
		file >> pos[1];
		file >> pos[2];
		hand_feature_pos.push_back(pos);
		hand_feature_idx.push_back(idx);
		hand_feature_weight.push_back(weight);
		x++;
		file >> count;
	}
	return;
}

void readFaceFeaturePos(string p_path, vector<Eigen::Vector3d> & face_feature_pos, vector<vector<int>> & face_feature_idx, vector<vector<double>> & face_feature_weight)
{
	ifstream file(p_path);
	int count;
	for (int i = 0; i < 54; i++)
	{
		file >> count;
		face_feature_idx[i].resize(1);
		face_feature_weight[i].resize(1);
		file >> face_feature_idx[i][0];
		Eigen::Vector3d content;
		file >> content(0);
		file >> content(1);
		file >> content(2);
		face_feature_pos[i] = content;
		face_feature_weight[i][0] = 1.0;
	}
	return;
}

Eigen::VectorXd readShapeCoefficient(string p_path)
{
	Eigen::VectorXd x(10);
	ifstream l(p_path);
	for (int i = 0; i < 10; i++)
	{
		l >> x(i);
	}
	l.close();
	return x;
}

Eigen::VectorXd readMotionCoefficient(string p_path)
{
	Eigen::VectorXd x(69);
	ifstream l(p_path);
	for (int i = 0; i < 69; i++)
	{
		l >> x(i);
	}
	l.close();
	return x;
}

void read_global(char* p_path, Eigen::Matrix3d & globalrotate, Eigen::Matrix3d & fromPctosmpl, Eigen::Matrix3d & fromsmpltoPc, double &s, Eigen::Vector3d & trans, Eigen::Vector3d & global_trans)
{
	ifstream file(p_path);
	file >> globalrotate(0, 0);
	file >> globalrotate(0, 1);
	file >> globalrotate(0, 2);
	file >> globalrotate(1, 0);
	file >> globalrotate(1, 1);
	file >> globalrotate(1, 2);
	file >> globalrotate(2, 0);
	file >> globalrotate(2, 1);
	file >> globalrotate(2, 2);
	file >> fromPctosmpl(0, 0);
	file >> fromPctosmpl(0, 1);
	file >> fromPctosmpl(0, 2);
	file >> fromPctosmpl(1, 0);
	file >> fromPctosmpl(1, 1);
	file >> fromPctosmpl(1, 2);
	file >> fromPctosmpl(2, 0);
	file >> fromPctosmpl(2, 1);
	file >> fromPctosmpl(2, 2);
	file >> fromsmpltoPc(0, 0);
	file >> fromsmpltoPc(0, 1);
	file >> fromsmpltoPc(0, 2);
	file >> fromsmpltoPc(1, 0);
	file >> fromsmpltoPc(1, 1);
	file >> fromsmpltoPc(1, 2);
	file >> fromsmpltoPc(2, 0);
	file >> fromsmpltoPc(2, 1);
	file >> fromsmpltoPc(2, 2);
	file >> s;
	file >> trans(0);
	file >> trans(1);
	file >> trans(2);
	file >> global_trans(0);
	file >> global_trans(1);
	file >> global_trans(2);
	file.close();
	return;
}
void readAdaptcoefficient(char * p_path, Eigen::VectorXd & p_coefficient)
{
	ifstream file(p_path);
	for (int i = 0; i < p_coefficient.size(); i++)
	{
		file >> p_coefficient(i);
	}
	file.close();
	return;
}

void readShapecoefficient(char * p_path, Eigen::VectorXd & p_coefficient)
{
	ifstream file(p_path);
	for (int i = 0; i < 10; i++)
	{
		file >> p_coefficient(i);
	}
	file.close();
	return;
}

void readPosecoefficient(char * p_path,  Eigen::VectorXd & p_coefficient)
{
	ifstream file(p_path);
	for (int i = 0; i < 400; i++)
	{
		file >> p_coefficient(i);
	}
	file.close();
	return;
}

void readcoefficient(char * p_path, Eigen::VectorXd & s_coefficient, Eigen::VectorXd & p_coefficient)
{
	ifstream file(p_path);
	for (int i = 0; i < 10; i++)
	{
		file >> s_coefficient(i);
	}
	for (int i = 0; i < 69; i++)
	{
		file >> p_coefficient(i);
	}
	file.close();
	return;
}

//void testHumanModel()
//{
//	char *meshpath = new char[100];
//	sprintf_s(meshpath, 100, "F:\\scan\\male_without_mouth\\0000.obj");
//	PGMesh *pg = new PGMesh();
//	OpenMesh::IO::read_mesh(*pg, meshpath);
//	Eigen::VectorXd y(69);
//	Eigen::VectorXd x(10);
//	Eigen::VectorXd facecoef(200);
//	Eigen::VectorXd gesturecoef(100);
//	Eigen::VectorXd poseparm(400);
//	Eigen::Matrix3d globalrotate, fromPctosmpl, fromsmpltoPc;
//	double s;
//	Eigen::Vector3d  trans, global_trans;
//	FaceBodyModel * fM = new FaceBodyModel(10, 200, pg);
//	Eigen::VectorXd out(fM->nV * 3);
//	sprintf_s(meshpath, 100, "F:\\scan\\male_expressions\\DE\\newC.txt");
//	fM->loadExpressionBases(meshpath, 200);
//	sprintf_s(meshpath, 100, "F:\\project\\splocs-master\\newC.txt");
//	fM->loadPoseBases(meshpath, 400);
//	sprintf_s(meshpath, 100, "F:\\scan\\male_hand\\DE\\newC.txt");
//	fM->loadgestureBases(meshpath, 100);
//	sprintf_s(meshpath, 100, "F:\\scan\\smpl\\male_shape\\mean.txt");
//	fM->loadMeanShape(meshpath);
//	sprintf_s(meshpath, 100, "F:\\scan\\smpl\\male_shape\\shapebasis.txt");
//	fM->loadShapeBases(meshpath);
//	Eigen::MatrixXd deformedV(fM->nV, 3);
//	fM->readBodyFeatureidx("body_feature_idx.txt");
//	int a[17] = { 12, 16, 18, 20, 1, 4, 7, 17, 19, 21, 2, 5, 8, 24, 25, 26, 27 };
//	fM->selectbodyfeature.resize(17);
//	copy(a, a + 17, fM->selectbodyfeature.begin());
//	fM->body_feature_idx1.resize(17);
//	fM->body_feature_weights1.resize(17);
//	for (int i = 0; i < 17; i++)
//	{
//		fM->body_feature_idx1[i] = fM->body_feature_idx[a[i]];
//		fM->body_feature_weights1[i] = fM->body_feature_weights[a[i]];
//	}
//	std::auto_ptr<mns::MyMeshHandler> p_meshHandler(new mns::MyMeshHandler(meshpath, fM->body_feature_idx1, fM->body_feature_weights1));
//	Eigen::VectorXd tmpV(fM->nV * 3);
//	for (int frame = 5757; frame < 5758; frame++)
//	{
//		for (int i_th = 2; i_th < 3; i_th++)
//		{
//			sprintf_s(meshpath, 100, "data\\%d\\%d\\coef.txt", frame, i_th);
//			readcoefficient(meshpath, x, y);
//			sprintf_s(meshpath, 100, "data\\%d\\%d\\global.txt", frame, i_th);
//			read_global(meshpath, globalrotate, fromPctosmpl, fromsmpltoPc, s, trans, global_trans);
//			sprintf_s(meshpath, 100, "data\\%d\\%d\\anchorpoint.txt", frame, i_th);
//			//fM->spCM->loadAnchorData(meshpath);
//			//fM->spCM->presolve();
//			sprintf_s(meshpath, 100, "data\\%d\\%d\\bodyconset.txt", frame, i_th);
//			fM->body_feature_pos1.resize(19);
//			readBodyFeaturePos(meshpath, fM->body_feature_pos1);
//			fM->body_feature_pos.resize(17);
//			fM->body_feature_pos[0] = fM->body_feature_pos1[0];
//			for (int i = 1; i < 17; i++)
//			{
//				fM->body_feature_pos[i] = fM->body_feature_pos1[i + 2];
//			}
//			sprintf_s(meshpath, 100, "data\\%d\\%d\\handconset.txt", frame, i_th);
//			fM->hand_feature_idx.resize(42);
//			fM->hand_feature_weights.resize(42);
//			fM->hand_feature_pos.resize(42);
//			readHandFeaturePos(meshpath, fM->hand_feature_pos, fM->hand_feature_idx, fM->hand_feature_weights);
//			sprintf_s(meshpath, 100, "data\\%d\\%d\\faceconset.txt", frame, i_th);
//			fM->face_feature_idx.resize(54);
//			fM->face_feature_pos.resize(54);
//			fM->face_feature_weights.resize(54);
//			readFaceFeaturePos(meshpath, fM->face_feature_pos, fM->face_feature_idx, fM->face_feature_weights);
//			int b[54] = { 9160, 3030, 152, 7173, 9811, 9249, 7258, 8782, 444, 4250, 5242, 1871, 1574, 6445, 4778, 5237, 932, 4336, 4407, 1160, 3098, 1099, 4613, 4578, 3294, 1127, 1050, 870, 4313, 4510, 6551, 3131, 3146, 1216, 2952, 1909, 1976, 6303, 5141, 5401, 5499, 5619, 2794, 1569, 7629, 7690, 4759, 1244, 4167, 6256, 1435, 890, 3063, 4303};
//			/*for (int i = 0; i < 54; i++)
//			{
//				fM->face_feature_idx[54 + i].resize(1);
//				fM->face_feature_idx[54 + i][0] = b[i];
//				fM->face_feature_idx[54 + i][0] = b[i];
//				fM->face_feature_weights[54 + i].resize(1);
//				fM->face_feature_weights[54 + i][0] = 1.0;
//			}*/
//
//
//			poseparm = fM->getPoseParm(x, y);
//			//cout << poseparm << endl;
//			Eigen::MatrixXd V = fM->f_x(x, poseparm);
//			Eigen::VectorXd currentFeaturepos = fM->getBodyPosfromCurrentV(V);
//			//优化全局旋转参数R,t
//			fM->findRigidTransform(currentFeaturepos, globalrotate, global_trans);
//
//			Eigen::MatrixXd sRV = globalrotate * (V.transpose());
//			sRV.colwise() += global_trans;
//			Eigen::MatrixXd sRV1 = sRV.transpose();
//			Eigen::MatrixXd sRV2 = s * (fromsmpltoPc * sRV);
//			sRV2.colwise() += trans;
//			Eigen::MatrixXd sRV3 = sRV2.transpose();
//			sprintf_s(meshpath, 100, "data\\%d\\%d\\initial_global.off", frame, i_th);
//			fM->spCM->outMesh1(sRV3, meshpath);
//			sprintf_s(meshpath, 100, "data\\%d\\%d\\initial_reconstruction.off", frame, i_th);
//			fM->spCM->outMesh1(sRV1, meshpath);
//			p_meshHandler->fitmode = 0;
//			
//			p_meshHandler->HandlerInit(meshpath, fM->body_feature_idx1, fM->body_feature_weights1);
//			p_meshHandler->p_DeformHandler->GetFeaturePos(fM->body_feature_pos);
//			p_meshHandler->p_DeformHandler->fitmode = 0;
//			//已知未知顶点位置，求DE
//			p_meshHandler->DeformationLD(selected_vid_vec, deformedV, deformedV);
//			sprintf_s(meshpath, 100, "data\\%d\\%d\\1_th_deformed.off", frame, i_th);
//			fM->spCM->outMesh1(deformedV, meshpath);
//			sRV2 = s * (fromsmpltoPc * deformedV.transpose());
//			sRV2.colwise() += trans;
//			sRV3 = sRV2.transpose();
//			sprintf_s(meshpath, 100, "data\\%d\\%d\\1_th_global_deformed.off", frame, i_th);
//			ofstream filefile1("before.txt");
//			filefile1 << poseparm << endl;
//			filefile1.close();
//			fM->spCM->outMesh1(sRV3, meshpath);
//			for (int i = 0; i < fM->nV; i++)
//			{
//				tmpV.segment(3 * i + 0, 3) = deformedV.row(i);
//			}
//			Eigen::VectorXd cm = fM->spCM->computeDiEdge(tmpV);
//			//for (int iter = 0; iter < 6; iter++)
//			//{
//			//	ofstream filefile2("after.txt");
//			//	fM->optimizebodyfromde(cm, poseparm);
//			//	filefile2 << poseparm << endl;
//			//	filefile2.close();
//			//	//V = fM->f_x(x, poseparm);
//			//	p_meshHandler->p_DeformHandler->UpdateWntVecWithBasis3(cm);
//
//			//	//Eigen::VectorXd currentFeaturepos = fM->getBodyPosfromCurrentV(V);
//			//	//优化全局旋转参数R,t
//			//	//fM->findRigidTransform(currentFeaturepos, globalrotate, global_trans);
//			//	//sRV = globalrotate * (V.transpose());
//			//	//sRV.colwise() += global_trans;
//			//	sRV1 = V;
//			//	sRV2 = s * (fromsmpltoPc * sRV1.transpose());
//			//	sRV2.colwise() += trans;
//			//	sRV3 = sRV2.transpose();
//			//	sprintf_s(meshpath, 100, "data\\%d\\%d\\body_%d_th_global_FBMM_iter.off", frame, i_th, iter);
//			//	fM->spCM->outMesh1(sRV3, meshpath);
//			//	sprintf_s(meshpath, 100, "data\\%d\\%d\\body_%d_th_fBMM_iter.off", frame, i_th, iter);
//			//	fM->spCM->outMesh1(V, meshpath);
//			//	p_meshHandler->HandlerInit(meshpath, fM->body_feature_idx1, fM->body_feature_weights1);
//			//	p_meshHandler->p_DeformHandler->GetFeaturePos(fM->body_feature_pos);
//			//	p_meshHandler->DeformationLD(selected_vid_vec, V, V);
//			//	sprintf_s(meshpath, 100, "data\\%d\\%d\\body_%d_th_deform_iter.off", frame, i_th, iter);
//			//	fM->spCM->outMesh1(V, meshpath);
//
//			//	for (int i = 0; i < fM->nV; i++)
//			//	{
//			//		tmpV.segment(3 * i + 0, 3) = V.row(i);
//			//	}
//			//	fM->spCM->setMesh(tmpV);
//			//	fM->spCM->calculateFacesNormals();
//			//	fM->spCM->calculateFacesFrame();
//			//	fM->spCM->computeDiEdge();
//			//	cm = fM->spCM->DiEdgeDataMatrix.row(0);
//			//}
//
//
//			p_meshHandler->fitmode = 1;
//			//开始优化手部姿态
//			sprintf_s(meshpath, 100, "data\\%d\\%d\\body_5_th_deform_iter.off", frame, i_th);
//
//			PGMesh * pgh = new PGMesh();
//			OpenMesh::IO::read_mesh(*pgh, meshpath);
//			fM->initFace(pgh);
//			p_meshHandler->HandlerInit(meshpath, fM->hand_feature_idx, fM->hand_feature_weights);
//			p_meshHandler->p_DeformHandler->GetFeaturePos(fM->hand_feature_pos);
//			p_meshHandler->p_DeformHandler->fitmode = 1;
//			p_meshHandler->p_DeformHandler->setFixed(fM->body_feature_idx1, fM->body_feature_weights1, fM->body_feature_pos);
//			p_meshHandler->DeformationLD(selected_vid_vec, deformedV, deformedV);
//			sprintf_s(meshpath, 100, "data\\%d\\%d\\hand_1_th_deformed.off", frame, i_th);
//			fM->spCM->outMesh1(deformedV, meshpath);
//			for (int i = 0; i < fM->nV; i++)
//			{
//				tmpV.segment(3 * i + 0, 3) = deformedV.row(i);
//			}
//			fM->spCM->setMesh(tmpV);
//			fM->spCM->calculateFacesNormals();
//			fM->spCM->calculateFacesFrame();
//			fM->spCM->computeDiEdge();
//			cm = fM->spCM->DiEdgeDataMatrix.row(0);
//			for (int iter = 0; iter < 8; iter++)
//			{
//				fM->optimizegesturefromde(cm, gesturecoef);
//				p_meshHandler->p_DeformHandler->UpdateWntVecWithBasis3(cm);
//				p_meshHandler->DeformationLD(selected_vid_vec, V, V);
//				sprintf_s(meshpath, 100, "data\\%d\\%d\\hand_%d_th_iter_FBMM.off", frame, i_th, iter);
//				fM->spCM->outMesh1(V, meshpath);
//				for (int i = 0; i < fM->nV; i++)
//				{
//					tmpV.segment(3 * i + 0, 3) = V.row(i);
//				}
//				sRV = V;
//				sRV2 = s * (fromsmpltoPc * sRV.transpose());
//				sRV2.colwise() += trans;
//				sRV3 = sRV2.transpose();
//				sprintf_s(meshpath, 100, "data\\%d\\%d\\hand_%d_th_global_iter_FBMM.off", frame, i_th, iter);
//				fM->spCM->outMesh1(sRV3, meshpath);
//				fM->spCM->setMesh(tmpV);
//				fM->spCM->calculateFacesNormals();
//				fM->spCM->calculateFacesFrame();
//				fM->spCM->computeDiEdge();
//				cm = fM->spCM->DiEdgeDataMatrix.row(0);
//			}
//
//			p_meshHandler->fitmode = 1;
//			//开始优化脸部表情
//			sprintf_s(meshpath, 100, "data\\%d\\%d\\body_5_th_deform_iter.off", frame, i_th);
//			
//			pgh= new PGMesh();
//			OpenMesh::IO::read_mesh(*pgh, meshpath);
//			fM->initFace(pgh);
//			p_meshHandler->HandlerInit(meshpath, fM->face_feature_idx, fM->face_feature_weights);
//			p_meshHandler->p_DeformHandler->GetFeaturePos(fM->face_feature_pos);
//			p_meshHandler->p_DeformHandler->fitmode = 1;
//			p_meshHandler->p_DeformHandler->setFixed(fM->body_feature_idx1, fM->body_feature_weights1, fM->body_feature_pos);
//			p_meshHandler->DeformationLD(selected_vid_vec, deformedV, deformedV);
//			sprintf_s(meshpath, 100, "data\\%d\\%d\\face_1_th_deformed.off", frame, i_th);
//			fM->spCM->outMesh1(deformedV, meshpath);
//			for (int i = 0; i < fM->nV; i++)
//			{
//				tmpV.segment(3 * i + 0, 3) = deformedV.row(i);
//			}
//			fM->spCM->setMesh(tmpV);
//			fM->spCM->calculateFacesNormals();
//			fM->spCM->calculateFacesFrame();
//			fM->spCM->computeDiEdge();
//			cm = fM->spCM->DiEdgeDataMatrix.row(0);
//			for (int iter = 0; iter < 5; iter++)
//			{
//				fM->optimizeexpressionfromde(cm, facecoef);
//				p_meshHandler->p_DeformHandler->UpdateWntVecWithBasis3(cm);
//				p_meshHandler->DeformationLD(selected_vid_vec, V, V);
//				sprintf_s(meshpath, 100, "data\\%d\\%d\\face_%d_th_iter_FBMM.off", frame, i_th, iter);
//				fM->spCM->outMesh1(V, meshpath);
//				for (int i = 0; i < fM->nV; i++)
//				{
//					tmpV.segment(3 * i + 0, 3) = V.row(i);
//				}
//				sRV = V;
//				sRV2 = s * (fromsmpltoPc * sRV.transpose());
//				sRV2.colwise() += trans;
//				sRV3 = sRV2.transpose();
//				sprintf_s(meshpath, 100, "data\\%d\\%d\\face_%d_th_global_iter_FBMM.off", frame, i_th, iter);
//				fM->spCM->outMesh1(sRV3, meshpath);
//				fM->spCM->setMesh(tmpV);
//				fM->spCM->calculateFacesNormals();
//				fM->spCM->calculateFacesFrame();
//				fM->spCM->computeDiEdge();
//				cm = fM->spCM->DiEdgeDataMatrix.row(0);
//			}
//			
//			
//
//			
//
//		}
//	}
//
//	//cout << y;
//}
void generateTwistingPose(int argc, char** argv)
{
	char *meshpath = new char[100];
	sprintf_s(meshpath, 100, "F:\\scan\\male_without_mouth\\0000.obj");
	PGMesh *pg = new PGMesh();
	OpenMesh::IO::read_mesh(*pg, meshpath);
	Eigen::VectorXd y(69);
	Eigen::VectorXd x(10);
	Eigen::VectorXd facecoef(200);
	Eigen::VectorXd gesturecoef(100);
	Eigen::VectorXd poseparm(400);
	Eigen::Matrix3d globalrotate, fromPctosmpl, fromsmpltoPc;
	double s;
	std::vector<int> layer = { 0, 1, 2, 3, 4, 5,6,7,8,9,10,11,12,13,14,15,16,17,18, 19, 20, 21, 22, 23 };
	std::vector<int> transvec = { 0, 6, 7, 8, 12, 13, 14 };
	std::vector<int> transvec1 = { 5, 8, 11, 14 };
	std::vector<int> layer1 = { 0, 1, 2, 3, 6, 9, 13, 14, 16, 17, 12, 18, 19, 4, 5 };
	std::vector<int> layer2 = { 18, 19, 4, 5 };
	std::vector<int> layer3 = { 7, 8, 10, 11, 20, 21, 22, 23 };
	Eigen::Vector3d  trans, global_trans;
	FaceBodyModel * fM = new FaceBodyModel(10, 200, pg);
	Eigen::VectorXd out(fM->nV * 3);
	sprintf_s(meshpath, 100, "F:\\scan\\male_expressions\\DE\\newC.txt");
	fM->loadExpressionBases(meshpath, 200);
	sprintf_s(meshpath, 100, "F:\\project\\splocs-master\\newC.txt");
	fM->loadPoseBases(meshpath, 400);
	sprintf_s(meshpath, 100, "F:\\scan\\male_hand\\DE\\newC.txt");
	fM->loadgestureBases(meshpath, 100);
	sprintf_s(meshpath, 100, "F:\\scan\\smpl\\male_shape\\mean.txt");
	fM->loadMeanShape(meshpath);
	sprintf_s(meshpath, 100, "F:\\scan\\smpl\\male_shape\\shapebasis.txt");
	fM->loadShapeBases(meshpath);
	Eigen::MatrixXd deformedV(fM->nV, 3);
	fM->readBodyFeatureidx("body_feature_idx.txt");
	int a[24] = { 0, 1, 2, 3, 4, 5, 6, 7,8 ,9,10, 11, 12, 13, 14, 15,16, 17, 18, 19, 20, 21, 22, 23 };
	fM->selectbodyfeature.resize(24);
	copy(a, a + 24, fM->selectbodyfeature.begin());
	fM->body_feature_idx1.resize(24);
	fM->body_feature_weights1.resize(24);
	for (int i = 0; i < 24; i++)
	{
		fM->body_feature_idx1[i] = fM->body_feature_idx[a[i]];
		fM->body_feature_weights1[i] = fM->body_feature_weights[a[i]];
	}
	std::auto_ptr<mns::MyMeshHandler> p_meshHandler(new mns::MyMeshHandler(meshpath, fM->body_feature_idx1, fM->body_feature_weights1));
	Eigen::VectorXd tmpV(fM->nV * 3);
	std::vector<int> smplbody0 = { 3500, 4985, 912, 5333, 1294, 3022, 2967, 6427, 3107, 6530,1619, 1665, 5212, 5165, 4515, 1029, 4486, 1001 };
	std::vector<int> fbmmbody0 = { 4785, 4956, 912, 5306, 1294, 2995, 2940, 6373, 3080, 6476 ,1619, 1665, 5185, 5138, 4488, 1029, 4590, 1001 };
	std::vector<int> smplbody1 = { 3500, 4985, 912, 5333, 1294, 3022, 2967, 6427, 3107, 6530 };
	std::vector<int> fbmmbody1 = { 4785, 4956, 912, 5306, 1294, 2995, 2940, 6373, 3080, 6476 };
	std::vector<int> smplbody2 = { 1619, 1665, 5212, 5165, 4515, 1029, 4486, 1001 };
	std::vector<int> fbmmbody2 = { 1619, 1665, 5185, 5138, 4488, 1029, 4590, 1001 };
	std::vector<int> smplbody3 = { 1970, 1930, 5669, 5560, 6723, 6581, 3323, 3203 };
	std::vector<int> fbmmbody3 = { 1970, 1930, 5642, 5533, 6669, 6527, 3296, 3176 };
	std::vector<int> smplidx(fM->nV);
	//ifstream smplidxfile("smpl_idx.txt");
	for (int i = 0; i < fM->nV; i++)
	{
		//double b;
		//smplidxfile >> b;
		smplidx[i] = i;
		//std::cout << smplidx[i] << endl;
	}
	//smplidxfile.close();
	std::vector<int> smpltag = {6899, 3721, 3743, 8707, 7077, 10475, 9665, 1326, 6418, 4717, 5338, 4854, 5402, 5409, 5388, 5386, 5394, 5437, 5657, 5614, 5528, 1408, 1379, 779, 1413, 1685, 1689, 1552, 1558, 1549, 2223, 2170, 2079, 2094, 1988, 2268};
	std::vector<int> fbmmtag(smpltag.size());
	std::vector<std::vector<int>> fbmmtagidx(smpltag.size());
	std::vector<std::vector<double>> fbmmtagweights(smpltag.size());
	std::vector<Eigen::Vector3d> fbmmtagpos(smpltag.size());
	std::vector<int> smpltag1 = {5701, 5556, 5814, 5623, 5622, 5689, 5530, 5503, 5658, 5431, 5438, 5618, 5619, 5823, 5501, 5799, 5951, 6023, 5628, 5609, 5563, 5714, 5891, 5793, 5949, 6066, 6076, 6172, 6151, 5722, 2221, 2000, 2001, 2198, 2135, 2205, 2312, 2422, 2385, 2542, 2495, 2617, 2636, 2106, 1938, 1960, 2698, 2724, 2287, 2080, 2092, 2076, 1983, 1988, 2127, 2057 ,2777, 2514, 2454, 2343, 2740, 2632};
	std::vector<int> fbmmtag1(smpltag1.size());
	std::vector<std::vector<int>> fbmmtagidx1(smpltag1.size());
	std::vector<std::vector<double>> fbmmtagweights1(smpltag1.size());
	std::vector<Eigen::Vector3d> fbmmtagpos1(smpltag1.size());
	for (int i = 0; i < smpltag.size(); i++)
	{
		fbmmtag[i] = smplidx[smpltag[i]];
		std::vector<int> idx = { fbmmtag[i] };
		std::vector<double> weight = { 1.0 };
		fbmmtagidx[i] = idx;
		fbmmtagweights[i] = weight;
	}
	for (int i = 0; i < smpltag1.size(); i++)
	{
		fbmmtag1[i] = smplidx[smpltag1[i]];
		std::vector<int> idx = { fbmmtag1[i] };
		std::vector<double> weight = { 1.0 };
		fbmmtagidx1[i] = idx;
		fbmmtagweights1[i] = weight;
	}
	PGMesh * tmpmesh = new PGMesh();
	for (int frame = 1111; frame < 1112; frame++)
	{

		sprintf_s(meshpath, 100, "twistingdata\\%d\\%04d.obj", frame, frame);
		OpenMesh::IO::read_mesh(*tmpmesh, meshpath);
		//OpenMesh::IO::write_mesh(*tmpmesh, "refMesh2.off");
		Eigen::MatrixXd smplposMatrix = fM->getVertexMatrix1(tmpmesh);
		delete tmpmesh;

		for (int i = 0; i < smpltag.size(); i++)
		{
			fbmmtagpos[i] = smplposMatrix.row(smpltag[i]).transpose();
		}
		for (int i = 0; i < smpltag1.size(); i++)
		{
			fbmmtagpos1[i] = smplposMatrix.row(smpltag1[i]).transpose();
		}
		sprintf_s(meshpath, 100, "twistingdata\\%d\\coef.txt", frame);
		readcoefficient(meshpath, x, y);
		poseparm = fM->getPoseParm(x, y);
		Eigen::VectorXd ref = fM->generateShapeWithoutPose(x);
		Eigen::MatrixXd refMatrix1(fM->nV, 3);
		for (int i = 0; i < fM->nV; i++)
		{
			refMatrix1.row(i) = ref.segment(3 * i, 3);
		}

		Eigen::MatrixXd V = fM->f_x(x, poseparm);
		sprintf_s(meshpath, 100, "twistingdata\\%d\\out1.off", frame);
		fM->spCM->outMesh1(V, meshpath);
		fM->body_feature_pos = fM->getJointsPos(x, y);
		Eigen::VectorXd currentFeaturepos = fM->getBodyPosfromCurrentV(V, layer);
		fM->findRigidTransform(currentFeaturepos, globalrotate, global_trans, fM->body_feature_pos);
		Eigen::MatrixXd refMatrix = (globalrotate * V.transpose()).transpose();
		refMatrix.rowwise() += global_trans.transpose();
		fM->spCM->outMesh1(refMatrix, "refMesh2.off");
		refMatrix1 = (globalrotate * refMatrix1.transpose()).transpose();
		refMatrix1.rowwise() += global_trans.transpose();
		fM->spCM->outMesh1(refMatrix1, "refmesh.off");

		ofstream outfile1("initialize.txt");
		for (int i = 0; i < 400; i++)
		{
			outfile1 << poseparm(i) << endl;
		}
		
		mns::MyMeshHandler* myH = new mns::MyMeshHandler("refmesh.off", "parameters5.txt", argc, argv);
		myH->p_DeformHandler->ReceiveBasis2("C.txt");
		std::vector<integer> v_id;
		std::vector<std::vector<int>> layer1_idx(layer1.size());
		std::vector<std::vector<double>> layer1_weights(layer1.size());
		std::vector<Eigen::Vector3d> layer1_pos(layer1.size());

		std::vector<std::vector<int>> layer2_idx(layer2.size());
		std::vector<std::vector<double>> layer2_weights(layer2.size());
		std::vector<Eigen::Vector3d> layer2_pos(layer2.size());

		std::vector<std::vector<int>> layer3_idx(layer3.size());
		std::vector<std::vector<double>> layer3_weights(layer3.size());
		std::vector<Eigen::Vector3d> layer3_pos(layer3.size());
		for (int i = 0; i < layer1.size(); i++)
		{
			std::vector<int> idx = fM->body_feature_idx[layer1[i]];
			std::vector<double> weight = fM->body_feature_weights[layer1[i]];
			Eigen::Vector3d pos = fM->body_feature_pos[layer1[i]];
			layer1_idx[i] = idx;
			layer1_weights[i] = weight;
			layer1_pos[i] = pos;
		}
		for (int i = 0; i < layer2.size(); i++)
		{
			std::vector<int> idx = fM->body_feature_idx[layer2[i]];
			std::vector<double> weight = fM->body_feature_weights[layer2[i]];
			Eigen::Vector3d pos = fM->body_feature_pos[layer2[i]];
			layer2_idx[i] = idx;
			layer2_weights[i] = weight;
			layer2_pos[i] = pos;
		}
		for (int i = 0; i < layer3.size(); i++)
		{
			std::vector<int> idx = fM->body_feature_idx[layer3[i]];
			std::vector<double> weight = fM->body_feature_weights[layer3[i]];
			Eigen::Vector3d pos = fM->body_feature_pos[layer3[i]];
			layer3_idx[i] = idx;
			layer3_weights[i] = weight;
			layer3_pos[i] = pos;
		}

		//第一层
		myH->p_DeformHandler->DeformInit(v_id, V, V);
		myH->p_DeformHandler->setconset(layer1_idx, layer1_weights, layer1_pos);
		myH->p_DeformHandler->DeformInit1(v_id, V, V);
		myH->p_DeformHandler->econ_threshold = 4.0;
		myH->p_DeformHandler->Deform(v_id, V, V);
		myH->p_DeformHandler->ModVMat(V);
		sprintf_s(meshpath, 100, "twistingdata\\%d\\layer1.off", frame);
		fM->spCM->outMesh1(V, meshpath);
		myH->p_DeformHandler->addconset(layer2_idx, layer2_weights, layer2_pos);
		myH->p_DeformHandler->DeformInit1(v_id, V, V);
		myH->p_DeformHandler->Deform(v_id, V, V);
		myH->p_DeformHandler->ModVMat(V);
		sprintf_s(meshpath, 100, "twistingdata\\%d\\layer2.off", frame);
		fM->spCM->outMesh1(V, meshpath);


		myH->p_DeformHandler->addconset(layer3_idx, layer3_weights, layer3_pos);
		myH->p_DeformHandler->DeformInit1(v_id, V, V);
		myH->p_DeformHandler->Deform(v_id, V, V);
		myH->p_DeformHandler->ModVMat(V);
		sprintf_s(meshpath, 100, "twistingdata\\%d\\layer3.off", frame);
		fM->spCM->outMesh1(V, meshpath);

		myH->p_DeformHandler->alpha_c = 8000;
		myH->p_DeformHandler->addconset(fbmmtagidx, fbmmtagweights, fbmmtagpos);
		myH->p_DeformHandler->DeformInit1(v_id, V, V);
		myH->p_DeformHandler->econ_threshold = 0.1;
		myH->p_DeformHandler->Deform(v_id, V, V);
		myH->p_DeformHandler->ModVMat(V);
		sprintf_s(meshpath, 100, "twistingdata\\%d\\layer4.off", frame);
		fM->spCM->outMesh1(V, meshpath);

		myH->p_DeformHandler->addconset(fbmmtagidx1, fbmmtagweights1, fbmmtagpos1);
		myH->p_DeformHandler->DeformInit1(v_id, V, V);
		myH->p_DeformHandler->econ_threshold = 0.1;
		myH->p_DeformHandler->Deform(v_id, V, V);
		myH->p_DeformHandler->ModVMat(V);
		sprintf_s(meshpath, 100, "twistingdata\\%d\\layer5.off", frame);
		fM->spCM->outMesh1(V, meshpath);
	}
}

void generateyuanqingtwistingpose(int argc, char** argv)
{
	char *meshpath = new char[100];
	sprintf_s(meshpath, 100, "F:\\scan\\male_without_mouth\\0000.obj");
	PGMesh *pg = new PGMesh();
	OpenMesh::IO::read_mesh(*pg, meshpath);
	Eigen::VectorXd y(69);
	Eigen::VectorXd x(10);
	Eigen::VectorXd facecoef(200);
	Eigen::VectorXd gesturecoef(100);
	Eigen::VectorXd poseparm(500);
	Eigen::Matrix3d globalrotate, fromPctosmpl, fromsmpltoPc;
	double s;
	std::vector<int> layer = { 0, 1, 2, 3, 4, 5,6,7,8,9,10,11,12,13,14,15,16,17,18, 19, 20, 21, 22, 23 };
	std::vector<int> transvec = { 0, 6, 7, 8, 12, 13, 14 };
	std::vector<int> transvec1 = { 5, 8, 11, 14 };
	std::vector<int> layer1 = { 0, 1, 2, 3, 6, 9, 13, 14, 16, 17, 12, 18, 19, 4, 5 };
	std::vector<int> layer2 = { 18, 19, 4, 5 };
	std::vector<int> layer3 = { 7, 8, 10, 11, 20, 21, 22, 23 };
	Eigen::Vector3d  trans, global_trans;
	FaceBodyModel * fM = new FaceBodyModel(10, 250, pg);
	Eigen::VectorXd out(fM->nV * 3);
	sprintf_s(meshpath, 100, "F:\\scan\\male_expressions\\DE\\newC.txt");
	fM->loadExpressionBases(meshpath, 200);
	/*sprintf_s(meshpath, 100, "F:\\project\\splocs-master\\newC.txt");
	fM->loadPoseBases(meshpath, 400);*/
	sprintf_s(meshpath, 100, "E:\\reconstruction\\250\\work\\C.txt");
	fM->loadPoseBases(meshpath, 500);
	sprintf_s(meshpath, 100, "F:\\scan\\male_hand\\DE\\newC.txt");
	fM->loadgestureBases(meshpath, 100);
	sprintf_s(meshpath, 100, "F:\\scan\\smpl\\male_shape\\mean.txt");
	fM->loadMeanShape(meshpath);
	sprintf_s(meshpath, 100, "F:\\scan\\smpl\\male_shape\\shapebasis.txt");
	fM->loadShapeBases(meshpath);
	Eigen::MatrixXd deformedV(fM->nV, 3);
	fM->readBodyFeatureidx("body_feature_idx.txt");
	int a[24] = { 0, 1, 2, 3, 4, 5, 6, 7,8 ,9,10, 11, 12, 13, 14, 15,16, 17, 18, 19, 20, 21, 22, 23 };
	fM->selectbodyfeature.resize(24);
	copy(a, a + 24, fM->selectbodyfeature.begin());
	fM->body_feature_idx1.resize(24);
	fM->body_feature_weights1.resize(24);
	for (int i = 0; i < 24; i++)
	{
		fM->body_feature_idx1[i] = fM->body_feature_idx[a[i]];
		fM->body_feature_weights1[i] = fM->body_feature_weights[a[i]];
	}
	sprintf_s(meshpath, 100, "F:\\scan\\male_without_mouth\\0000.obj");
	std::auto_ptr<mns::MyMeshHandler> p_meshHandler(new mns::MyMeshHandler(meshpath, fM->body_feature_idx1, fM->body_feature_weights1));
	Eigen::VectorXd tmpV(fM->nV * 3);
	std::vector<int> smplbody0 = { 3500, 4985, 912, 5333, 1294, 3022, 2967, 6427, 3107, 6530,1619, 1665, 5212, 5165, 4515, 1029, 4486, 1001 };
	std::vector<int> fbmmbody0 = { 4785, 4956, 912, 5306, 1294, 2995, 2940, 6373, 3080, 6476 ,1619, 1665, 5185, 5138, 4488, 1029, 4590, 1001 };
	std::vector<int> smplbody1 = { 3500, 4985, 912, 5333, 1294, 3022, 2967, 6427, 3107, 6530 };
	std::vector<int> fbmmbody1 = { 4785, 4956, 912, 5306, 1294, 2995, 2940, 6373, 3080, 6476 };
	std::vector<int> smplbody2 = { 1619, 1665, 5212, 5165, 4515, 1029, 4486, 1001 };
	std::vector<int> fbmmbody2 = { 1619, 1665, 5185, 5138, 4488, 1029, 4590, 1001 };
	std::vector<int> smplbody3 = { 1970, 1930, 5669, 5560, 6723, 6581, 3323, 3203 };
	std::vector<int> fbmmbody3 = { 1970, 1930, 5642, 5533, 6669, 6527, 3296, 3176 };
	std::vector<int> smplidx(fM->nV);
	//ifstream smplidxfile("smpl_idx.txt");
	for (int i = 0; i < fM->nV; i++)
	{
		//double b;
		//smplidxfile >> b;
		smplidx[i] = i;
		//std::cout << smplidx[i] << endl;
	}
	//smplidxfile.close();
	std::vector<int> smpltag = { 6899, 3721, 3743, 8707, 7077, 10475, 9665, 1326, 6418, 4717, 5338, 4854, 5402, 5409, 5388, 5386, 5394, 5437, 5657, 5614, 5528, 1408, 1379, 779, 1413, 1685, 1689, 1552, 1558, 1549, 2223, 2170, 2079, 2094, 1988, 2268 };
	std::vector<int> fbmmtag(smpltag.size());
	std::vector<std::vector<int>> fbmmtagidx(smpltag.size());
	std::vector<std::vector<double>> fbmmtagweights(smpltag.size());
	std::vector<Eigen::Vector3d> fbmmtagpos(smpltag.size());
	std::vector<int> smpltag1 = { 5701, 5556, 5814, 5623, 5622, 5689, 5530, 5503, 5658, 5431, 5438, 5618, 5619, 5823, 5501, 5799, 5951, 6023, 5628, 5609, 5563, 5714, 5891, 5793, 5949, 6066, 6076, 6172, 6151, 5722, 2221, 2000, 2001, 2198, 2135, 2205, 2312, 2422, 2385, 2542, 2495, 2617, 2636, 2106, 1938, 1960, 2698, 2724, 2287, 2080, 2092, 2076, 1983, 1988, 2127, 2057 ,2777, 2514, 2454, 2343, 2740, 2632 };
	std::vector<int> fbmmtag1(smpltag1.size());
	std::vector<std::vector<int>> fbmmtagidx1(smpltag1.size());
	std::vector<std::vector<double>> fbmmtagweights1(smpltag1.size());
	std::vector<Eigen::Vector3d> fbmmtagpos1(smpltag1.size());
	for (int i = 0; i < smpltag.size(); i++)
	{
		fbmmtag[i] = smplidx[smpltag[i]];
		std::vector<int> idx = { fbmmtag[i] };
		std::vector<double> weight = { 1.0 };
		fbmmtagidx[i] = idx;
		fbmmtagweights[i] = weight;
	}
	for (int i = 0; i < smpltag1.size(); i++)
	{
		fbmmtag1[i] = smplidx[smpltag1[i]];
		std::vector<int> idx = { fbmmtag1[i] };
		std::vector<double> weight = { 1.0 };
		fbmmtagidx1[i] = idx;
		fbmmtagweights1[i] = weight;
	}
	PGMesh * tmpmesh = new PGMesh();
	for (int frame = 1111; frame < 1112; frame++)
	{
		std::vector<std::vector<int>> layer1_idx;
		std::vector<std::vector<double>> layer1_weights;
		std::vector<Eigen::Vector3d> layer1_pos;
		sprintf_s(meshpath, 100, "twistingdata\\%d\\body_conset.txt", frame);
		readHandFeaturePos(meshpath, layer1_pos, layer1_idx, layer1_weights);
		sprintf_s(meshpath, 100, "twistingdata\\%d\\global.txt", frame);
		read_global(meshpath, globalrotate, fromPctosmpl, fromsmpltoPc, s, trans, global_trans);
		std::vector<std::vector<int>> layer2_idx;
		std::vector<std::vector<double>> layer2_weights;
		std::vector<Eigen::Vector3d> layer2_pos;
		sprintf_s(meshpath, 100, "twistingdata\\%d\\lower_conset.txt", frame);
		readHandFeaturePos(meshpath, layer2_pos, layer2_idx, layer2_weights);

		std::vector<std::vector<int>> layer3_idx;
		std::vector<std::vector<double>> layer3_weights;
		std::vector<Eigen::Vector3d> layer3_pos;
		sprintf_s(meshpath, 100, "twistingdata\\%d\\upper_conset.txt", frame);
		readHandFeaturePos(meshpath, layer3_pos, layer3_idx, layer3_weights);
		//OpenMesh::IO::write_mesh(*tmpmesh, "refMesh2.off");
		
		sprintf_s(meshpath, 100, "twistingdata\\%d\\coef.txt", frame);
		readcoefficient(meshpath, x, y);
		poseparm = fM->getPoseParm(x, y);
		Eigen::VectorXd ref = fM->generateShapeWithoutPose(x);
		Eigen::MatrixXd refMatrix1(fM->nV, 3);
		for (int i = 0; i < fM->nV; i++)
		{
			refMatrix1.row(i) = ref.segment(3 * i, 3);
		}

		Eigen::MatrixXd V = fM->f_x(x, poseparm);
		sprintf_s(meshpath, 100, "twistingdata\\%d\\out1.off", frame);
		fM->spCM->outMesh1(V, meshpath);
		fM->body_feature_pos = fM->getJointsPos(x, y);
		Eigen::VectorXd currentFeaturepos = fM->getBodyPosfromCurrentV(V, layer1_idx, layer1_weights);
		fM->findRigidTransform(currentFeaturepos, globalrotate, global_trans, layer1_pos);
		Eigen::MatrixXd refMatrix = (globalrotate * V.transpose()).transpose();
		refMatrix.rowwise() += global_trans.transpose();
		//fM->spCM->outMesh1(refMatrix, "refMesh2.off");
		refMatrix1 = (globalrotate * refMatrix1.transpose()).transpose();
		refMatrix1.rowwise() += global_trans.transpose();
		fM->spCM->outMesh1(refMatrix, "refmesh.off");
		fM->spCM->outMesh1(refMatrix1, "refMesh2.off");
		ofstream outfile1("initialize.txt");
		for (int i = 0; i < 400; i++)
		{
			outfile1 << poseparm(i) << endl;
		}
		sprintf_s(meshpath, 100, "twistingdata\\%d\\handconset.txt", frame);
		readHandFeaturePos(meshpath, fM->hand_feature_pos, fM->hand_feature_idx, fM->hand_feature_weights);
		for (int h_idx = 0; h_idx < fM->hand_feature_pos.size(); h_idx++)
		{
			Eigen::Vector3d v0 = fM->hand_feature_pos[h_idx];
			Eigen::Vector3d v1 = (1 / s) * (fromPctosmpl * (v0 - trans));
			fM->hand_feature_pos[h_idx] = v1;
			cout << v1(0) << " " << v1(1) << " " << v1(2) << endl;
		}
		mns::MyMeshHandler* myH = new mns::MyMeshHandler("refmesh.off", "parameters1.txt", argc, argv);
		//myH->p_DeformHandler->ReceiveBasis2("C.txt");
		myH->p_DeformHandler->ReceiveBasis2("E:\\reconstruction\\250\\work\\newC.txt");
		std::vector<integer> v_id;
		
		

		//第一层
		//myH->p_DeformHandler->DeformInit(v_id, V, V);
		//myH->p_DeformHandler->setconset(layer1_idx, layer1_weights, layer1_pos);
		//myH->p_DeformHandler->addconset(layer3_idx, layer3_weights, layer3_pos);
		//myH->p_DeformHandler->addconset(layer2_idx, layer2_weights, layer2_pos);
		//
		//myH->p_DeformHandler->DeformInit1(v_id, V, V);
		//myH->p_DeformHandler->econ_threshold = 4.0;
		//myH->p_DeformHandler->Deform(v_id, V, V);
		//myH->p_DeformHandler->ModVMat(V);
		//sprintf_s(meshpath, 100, "twistingdata\\%d\\layer1.off", frame);
		//fM->spCM->outMesh1(V, meshpath);
		//fM->spCM->outMesh1(V, "refmesh2.off");
		//myH = new mns::MyMeshHandler("refmesh2.off", "parameters5.txt", argc, argv);
		////myH->p_DeformHandler->ReceiveBasis2("C.txt");
		//myH->p_DeformHandler->ReceiveBasis2("E:\\reconstruction\\250\\work\\newC.txt");
		//myH->p_DeformHandler->DeformInit(v_id, V, V);
		//myH->p_DeformHandler->setconset(layer1_idx, layer1_weights, layer1_pos);
		//myH->p_DeformHandler->addconset(layer2_idx, layer2_weights, layer2_pos);
		//myH->p_DeformHandler->addconset(layer3_idx, layer3_weights, layer3_pos);
		//myH->p_DeformHandler->DeformInit1(v_id, V, V);
		//myH->p_DeformHandler->Deform(v_id, V, V);
		//myH->p_DeformHandler->ModVMat(V);
		//sprintf_s(meshpath, 100, "twistingdata\\%d\\layer2.off", frame);
		//fM->spCM->outMesh1(V, meshpath);

		//fM->spCM->outMesh1(V, "refmesh2.off");
		//myH = new mns::MyMeshHandler("refmesh2.off", "parameters5.txt", argc, argv);
		////myH->p_DeformHandler->ReceiveBasis2("C.txt");
		//myH->p_DeformHandler->ReceiveBasis2("E:\\reconstruction\\250\\work\\newC.txt");
		//myH->p_DeformHandler->DeformInit(v_id, V, V);
		//myH->p_DeformHandler->alpha_c = 2000;
		//myH->p_DeformHandler->econ_threshold = 20.0;
		//myH->p_DeformHandler->setconset(layer1_idx, layer1_weights, layer1_pos);
		//myH->p_DeformHandler->addconset(layer2_idx, layer2_weights, layer2_pos);
		//myH->p_DeformHandler->addconset(layer3_idx, layer3_weights, layer3_pos);
		//myH->p_DeformHandler->DeformInit1(v_id, V, V);
		//myH->p_DeformHandler->Deform(v_id, V, V);
		//myH->p_DeformHandler->ModVMat(V);
		//sprintf_s(meshpath, 100, "twistingdata\\%d\\layer3.off", frame);
		//fM->spCM->outMesh1(V, meshpath);

		//myH->p_DeformHandler->alpha_c = 8000;
		//myH->p_DeformHandler->addconset(fM->hand_feature_idx, fM->hand_feature_weights, fM->hand_feature_pos);
		//myH->p_DeformHandler->DeformInit1(v_id, V, V);
		//myH->p_DeformHandler->econ_threshold = 0.1;
		//myH->p_DeformHandler->Deform(v_id, V, V);
		//myH->p_DeformHandler->ModVMat(V);
		//sprintf_s(meshpath, 100, "twistingdata\\%d\\layer4.off", frame);
		//fM->spCM->outMesh1(V, meshpath);

		//手部姿态
		PGMesh * pg1 = new PGMesh();
		outfile1 = ofstream("initialize.txt");
		for (int i = 0; i < 100; i++)
		{
			outfile1 << 0 << endl;
		}
		outfile1.close();
		sprintf_s(meshpath, 100, "twistingdata\\%d\\layer4.off", frame);
		OpenMesh::IO::read_mesh(*pg1, meshpath);
		Eigen::VectorXd tempV1 = fM->getVertexMatrix(pg1);
		delete pg1;
		Eigen::MatrixXd V10(3, fM->nV);
		for (int i = 0; i < fM->nV; i++)
		{
			V10.col(i) = tempV1.segment(3 * i, 3);
		}
		
		Eigen::MatrixXd V11 = ( V10).transpose();
		//sRV = sRV2.transpose();

		fM->spCM->outMesh1(V11, "refMesh2.off");
		Eigen::MatrixXd V12(fM->nV, 3);
		myH = new mns::MyMeshHandler("refMesh2.off", "parameters3.txt", argc, argv);
		myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\male_hand\\DE\\newC.txt");
		myH->p_DeformHandler->DeformInit(v_id, V12, V12);
		std::vector<int> handfix = { 4087, 5177, 3038, 2896, 3523, 1678, 1668, 10450, 5221, 4840, 1709 };
		std::vector<std::vector<int>> handfixidx;
		std::vector<std::vector<double>> handfixweights;
		std::vector<Eigen::Vector3d> handfixpos;
		for (int i = 0; i < handfix.size(); i++)
		{
			vector<int> idx = { handfix[i] };
			vector<double> weight = { 1.0 };
			Eigen::Vector3d coord;
			coord = V11.row(handfix[i]).transpose();
			handfixidx.push_back(idx);
			handfixpos.push_back(coord);
			handfixweights.push_back(weight);
		}
		myH->p_DeformHandler->setconset(handfixidx, handfixweights, handfixpos);
		/*myH->p_DeformHandler->setconset(layer1_idx, layer1_weights, layer1_pos);
		myH->p_DeformHandler->addconset(layer2_idx, layer2_weights, layer2_pos);
		myH->p_DeformHandler->addconset(layer3_idx, layer3_weights, layer3_pos);*/
		myH->p_DeformHandler->addconset(fM->hand_feature_idx, fM->hand_feature_weights, fM->hand_feature_pos);

		//myH->p_DeformHandler->DeformInit1(v_id, V, V);
		
		//if (newstate == 1)
		{
			myH->p_DeformHandler->Deform(v_id, V12, V12);
		}

		
		myH->p_DeformHandler->ModVMat(V12);

		Eigen::MatrixXd V13 = s * (fromsmpltoPc * V12.transpose());
		V13.colwise() += trans;
		Eigen::MatrixXd V14 = V13.transpose();
		//Eigen::MatrixXd V14 = V12;
		sprintf_s(meshpath, 100, "twistingdata\\%d\\layer5.off", frame);
		fM->spCM->outMesh1(V14, meshpath);
	}
}

void generateDatafromModel(int argc, char** argv)
{
	std::cout << "start" << std::endl;
	char *meshpath = new char[100];
	sprintf_s(meshpath, 100, "F:\\scan\\male_without_mouth\\0000.obj");
	PGMesh *pg = new PGMesh();
	OpenMesh::IO::read_mesh(*pg, meshpath);
	Eigen::VectorXd y(69);//姿态(69个参数)
	Eigen::VectorXd x(10);//10个形状系数 
	x.setZero();
	Eigen::VectorXd facecoef(200);
	Eigen::VectorXd gesturecoef(100);
	Eigen::VectorXd poseparm(400);
	FaceBodyModel * fM = new FaceBodyModel(10, 200, pg);
	Eigen::VectorXd out(fM->nV * 3);
	sprintf_s(meshpath, 100, "F:\\scan\\male_expressions\\DE\\newC.txt");
	fM->loadExpressionBases(meshpath, 200);
	sprintf_s(meshpath, 100, "E:\\Graphics\\bachelor_thesis\\Code\\splocs-master\\newC.txt");//F:\\project\\splocs-master\\newC.txt
	fM->loadPoseBases(meshpath, 400);//no data
	sprintf_s(meshpath, 100, "F:\\scan\\male_hand\\DE\\newC.txt");
	fM->loadgestureBases(meshpath, 100);
	sprintf_s(meshpath, 100, "F:\\scan\\smpl\\male_shape\\mean.txt");
	fM->loadMeanShape(meshpath);
	sprintf_s(meshpath, 100, "F:\\scan\\smpl\\male_shape\\shapebasis.txt");
	fM->loadShapeBases(meshpath); 
	Eigen::MatrixXd deformedV(fM->nV, 3);
	fM->readBodyFeatureidx("body_feature_idx.txt");
	int a[19] = { 12, 15, 0, 16, 18, 20, 1, 4, 7, 17, 19, 21, 2, 5, 8, 24, 25, 26, 27 };
	fM->selectbodyfeature.resize(19);
	copy(a, a + 19, fM->selectbodyfeature.begin());
	fM->body_feature_idx1.resize(19);
	fM->body_feature_weights1.resize(19);
	for (int i = 0; i < 19; i++)
	{
		fM->body_feature_idx1[i] = fM->body_feature_idx[a[i]];
		fM->body_feature_weights1[i] = fM->body_feature_weights[a[i]];
	}
	for (int frameno = 1; frameno < 101; frameno++)
	{
		//sprintf_s(meshpath, 100, "generatedata\\%d\\posecoef.txt", frameno);
		//poseparm = fM->getPoseParm(x, y);
		//readPosecoefficient(meshpath, poseparm);
		//sprintf_s(meshpath, 100, "generatedata\\%d\\shapecoef.txt", frameno);
		//readShapecoefficient(meshpath, x);
		//sprintf_s(meshpath, 100, "generatedata\\%d\\posecoef.txt", frameno);
		//poseparm = fM->getPoseParm(x, y);
		sprintf_s(meshpath, 100, "E:\\Graphics\\bachelor_thesis\\Code\\smplify-x-master\\output_all_smpl\\results\\%02d_img\\coef.txt", frameno);
		readcoefficient(meshpath, x, y);
		poseparm = fM->getPoseParm(x, y);
		//std::cout << "Here is the vector x:\n" << x << std::endl;
		//std::cout << "Here is the vector y:\n" << poseparm << std::endl;
		//sprintf_s(meshpath, 100, "generatedata\\%d\\facecoef.txt", frameno);
		//readAdaptcoefficient(meshpath, facecoef);
		//sprintf_s(meshpath, 100, "generatedata\\%d\\handcoef.txt", frameno);
		//readAdaptcoefficient(meshpath, gesturecoef);
		//gesturecoef = gesturecoef * 2;
		facecoef.setZero();
		gesturecoef.setZero();
		Eigen::MatrixXd V = fM->f_x(x, poseparm, facecoef, gesturecoef);//**important**//
		sprintf_s(meshpath, 100, "generatedata\\EHF\\%02d_mesh_before.obj", frameno, frameno);//generatedata\\%d\\mesh1_%04d.obj
		fM->spCM->outMesh(V, meshpath);
		std::cout << frameno << std::endl;
	}
	std::cout << "end" << std::endl;
}

void testFaceRegistritation(int argc, char** argv)
{
	char *meshpath = new char[100];
	sprintf_s(meshpath, 100, "F:\\scan\\male_without_mouth\\0000.obj");
	PGMesh *pg = new PGMesh();
	OpenMesh::IO::read_mesh(*pg, meshpath);
	Eigen::VectorXd y(69);
	Eigen::VectorXd x(10);
	Eigen::VectorXd facecoef(200);
	Eigen::VectorXd gesturecoef(100);
	Eigen::VectorXd poseparm(400);
	Eigen::Matrix3d globalrotate, fromPctosmpl, fromsmpltoPc;
	Eigen::Vector3d global_trans;
	FaceBodyModel * fM = new FaceBodyModel(10, 200, pg);
	sprintf_s(meshpath, 100, "F:\\scan\\male_expressions\\DE\\newC.txt");
	fM->loadExpressionBases(meshpath, 200);
	fM->loadDenseCorrespondenceModule();
	sprintf_s(meshpath, 100, "F:\\project\\splocs-master\\newC.txt");
	fM->loadPoseBases(meshpath, 400);
	sprintf_s(meshpath, 100, "F:\\scan\\male_hand\\DE\\newC.txt");
	fM->loadgestureBases(meshpath, 100);
	sprintf_s(meshpath, 100, "F:\\scan\\smpl\\male_shape\\mean.txt");
	fM->loadMeanShape(meshpath);
	sprintf_s(meshpath, 100, "F:\\scan\\smpl\\male_shape\\shapebasis.txt");
	fM->loadShapeBases(meshpath);
	Eigen::MatrixXd deformedV(fM->nV, 3);
	fM->readBodyFeatureidx("body_feature_idx.txt");
	int a[19] = { 12, 15, 0, 16, 18, 20, 1, 4, 7, 17, 19, 21, 2, 5, 8, 24, 25, 26, 27 };
	fM->selectbodyfeature.resize(19);
	copy(a, a + 19, fM->selectbodyfeature.begin());
	fM->body_feature_idx1.resize(19);
	fM->body_feature_weights1.resize(19);
	for (int i = 0; i < 19; i++)
	{
		fM->body_feature_idx1[i] = fM->body_feature_idx[a[i]];
		fM->body_feature_weights1[i] = fM->body_feature_weights[a[i]];
	}
	PGMesh * tmpmesh = new PGMesh();
	for (int frame = 29; frame < 30; frame++)
	{
		fM->load_pointcloud(frame);
		for (int i_th = 14; i_th < 15; i_th++)
		{
			sprintf_s(meshpath, 100, "face_scan\\%d\\%d\\shape_coef.txt", frame, i_th);
			readcoefficient(meshpath, x, y);
			poseparm = fM->getPoseParm(x, y);

			//cout << poseparm << endl;
			Eigen::MatrixXd refMatrix = fM->f_x(x, poseparm);
			
			//优化全局旋转参数R,t
			//fM->findRigidTransform(currentFeaturepos, globalrotate, global_trans, fM->body_feature_pos1);
			
			//Eigen::VectorXd ref = fM->generateShapeWithoutPose(x);
			
			/*sprintf_s(meshpath, 100, "face_scan\\%d\\%d\\layer9.off", frame, i_th);
			OpenMesh::IO::read_mesh(*tmpmesh, meshpath);
			Eigen::MatrixXd refMatrix = fM->getVertexMatrix1(tmpmesh);*/
			
			fM->face_feature_idx.resize(0);
			fM->face_feature_pos.resize(0);
			fM->face_feature_weights.resize(0);
			sprintf_s(meshpath, 100, "face_scan\\%d\\%d\\face_conset.txt", frame, i_th);
			readFacecontourFeaturePos(meshpath, fM->face_feature_pos, fM->face_feature_idx, fM->face_feature_weights);
			Eigen::VectorXd currentFeaturepos = fM->getFacePosfromCurrentV(refMatrix);
			fM->findRigidTransform(currentFeaturepos, globalrotate, global_trans, fM->face_feature_pos);
			refMatrix = (globalrotate * refMatrix.transpose()).transpose();
			refMatrix.rowwise() += global_trans.transpose();
			fM->spCM->outMesh1(refMatrix, "refMesh2.off");
			sprintf_s(meshpath, 100, "face_scan\\%d\\%d\\layer3.off", frame, i_th);
			
			sprintf_s(meshpath, 100, "face_scan\\%d\\%d\\initial_reconstruction.off", frame, i_th);
			fM->spCM->outMesh1(refMatrix, meshpath);
			fM->spCM->outMesh1(refMatrix, "refMesh2.off");
			std::vector<integer> v_id;
			mns::MyMeshHandler* myH = new mns::MyMeshHandler("refMesh2.off", "parameters2.txt", argc, argv);
			myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\male_expressions\\DE\\C.txt");
			Eigen::MatrixXd V7(fM->nV, 3);
			myH->p_DeformHandler->DeformInit(v_id, V7, V7);
			std::vector<int> facefix = { 5298, 1864, 6448, 8387, 4782, 5105, 1651 };
			std::vector<std::vector<int>> facefixidx;
			std::vector<std::vector<double>> facefixweights;
			std::vector<Eigen::Vector3d> facefixpos;
			for (int i = 0; i < facefix.size(); i++)
			{
				vector<int> idx = { facefix[i] };
				vector<double> weight = { 1.0 };
				Eigen::Vector3d coord;
				coord = refMatrix.row(facefix[i]).transpose();
				facefixidx.push_back(idx);
				facefixpos.push_back(coord);
				facefixweights.push_back(weight);
			}
			myH->p_DeformHandler->setconset(facefixidx, facefixweights, facefixpos);
			myH->p_DeformHandler->addconset(fM->face_feature_idx, fM->face_feature_weights, fM->face_feature_pos);
			myH->p_DeformHandler->DeformInit1(v_id, V7, V7);
			myH->p_DeformHandler->econ_threshold = 1.0;
			myH->p_DeformHandler->Deform(v_id, V7, V7);
			myH->p_DeformHandler->ModVMat(V7);
			sprintf_s(meshpath, 100, "face_scan\\%d\\%d\\keypoint.off", frame, i_th);
			fM->spCM->outMesh1(V7, meshpath);
			sprintf_s(meshpath, 100, "face_scan\\%d\\%d\\layer9.off", frame, i_th);
			fM->spCM->outMesh1(V7, meshpath);
			sprintf_s(meshpath, 100, "face_scan\\%d\\%d\\layer9_0.off", frame, i_th);
			fM->spCM->outMesh1(V7, meshpath);
			std::vector<std::vector<int>> denceidx;
			std::vector<std::vector<double>> denceweights;
			std::vector<Eigen::Vector3d> dencepos;
			PGMesh *pg1;
			for (int lll = 0; lll < 4; lll++)
			{
				fM->findCoressPairs(frame, i_th, 9);
				sprintf_s(meshpath, 100, "face_scan\\%d\\%d\\coress.txt", frame, i_th);
				readCorespondence(meshpath, dencepos, denceidx, denceweights);
				pg1 = new PGMesh();
				sprintf_s(meshpath, 100, "face_scan\\%d\\%d\\layer9_%d.off", frame, i_th, lll);
				OpenMesh::IO::read_mesh(*pg1, meshpath);
				Eigen::VectorXd tempV1 = fM->getVertexMatrix(pg1);
				delete pg1;
				Eigen::MatrixXd V10(3, fM->nV);
				for (int i = 0; i < fM->nV; i++)
				{
					V10.col(i) = tempV1.segment(3 * i, 3);
				}
				
				Eigen::MatrixXd V12 = (V10).transpose();



				fM->spCM->outMesh1(V12, "refMesh2.off");
				sprintf_s(meshpath, 100, "refMesh2.off", frame, i_th);
				mns::MyMeshHandler *myH = new mns::MyMeshHandler(meshpath, "parameters2.txt", argc, argv);
				myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\male_expressions\\DE\\C.txt");

				myH->p_DeformHandler->DeformInit(v_id, V12, V12);

				//myH->p_DeformHandler->setconset(handfixidx, handfixweights, handfixpos);
				myH->p_DeformHandler->setconset(facefixidx, facefixweights, facefixpos);
				myH->p_DeformHandler->addconset(fM->face_feature_idx, fM->face_feature_weights, fM->face_feature_pos);
				myH->p_DeformHandler->addconset(denceidx, denceweights, dencepos);

				//myH->p_DeformHandler->DeformInit1(v_id, V, V);
				myH->p_DeformHandler->Deform(v_id, V12, V12);
				myH->p_DeformHandler->ModVMat(V12);



				sprintf_s(meshpath, 100, "face_scan\\%d\\%d\\layer9_%d.off", frame, i_th, lll+1);
				fM->spCM->outMesh1(V12, meshpath);
			}
		}
	}
}

void testHumanModel7(int argc, char** argv)
{
	char *meshpath = new char[100];
	sprintf_s(meshpath, 100, "F:\\scan\\male_without_mouth\\0000.obj");
	PGMesh *pg = new PGMesh();
	OpenMesh::IO::read_mesh(*pg, meshpath);
	Eigen::VectorXd y(69);
	Eigen::VectorXd x(10);
	Eigen::VectorXd facecoef(200);
	Eigen::VectorXd gesturecoef(100);
	Eigen::VectorXd poseparm(400);
	Eigen::Matrix3d globalrotate, fromPctosmpl, fromsmpltoPc;
	double s;
	std::vector<int> layer = { 0, 1, 2, 3, 4, 5,6,7,8,9,10,11,12,13,14,15,16,17,18,19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32 };
	std::vector<int> transvec = { 0, 6, 7, 8, 12, 13, 14 };
	std::vector<int> transvec1 = { 5, 8, 11, 14 };
	std::vector<int> layer1 = { 12, 13, 14, 16, 17, 6, 0, 1, 2, 4, 5, 18, 19, 24, 25, 26, 27, 28, 29, 30, 31, 32 };
	std::vector<int> layer2 = { 20, 21, 22, 23, 7, 8, 10, 11 };
	std::vector<int> layer3 = { 5, 8, 11, 14 };
	Eigen::Vector3d  trans, global_trans;
	FaceBodyModel * fM = new FaceBodyModel(10, 200, pg);
	//FaceBodyModel * fM = new FaceBodyModel(10, 250, pg);
	Eigen::VectorXd out(fM->nV * 3);
	sprintf_s(meshpath, 100, "F:\\scan\\male_expressions\\DE\\newC.txt");
	fM->loadExpressionBases(meshpath, 200);
	/*sprintf_s(meshpath, 100, "F:\\project\\splocs-master\\newC.txt");
	fM->loadPoseBases(meshpath, 400);*/
	/*sprintf_s(meshpath, 100, "F:\\scan\\MFaust\\scripts\\male\\mixedDE\\newC.txt");
	fM->loadPoseBases(meshpath, 500);*/
	sprintf_s(meshpath, 100, "F:\\scan\\male_without_mouth\\DE\\C.txt");//female->male
	fM->loadPoseBases(meshpath, 400);
	sprintf_s(meshpath, 100, "F:\\scan\\female_hand\\DE\\newC.txt");
	fM->loadgestureBases(meshpath, 100);
	sprintf_s(meshpath, 100, "F:\\scan\\smpl\\male_shape\\mean.txt");//female->male
	fM->loadMeanShape(meshpath);
	sprintf_s(meshpath, 100, "F:\\scan\\smpl\\male_shape\\shapebasis.txt");//female->male
	fM->loadShapeBases(meshpath);

	Eigen::MatrixXd deformedV(fM->nV, 3);
	fM->readBodyFeatureidx1("body_feature_idx1.txt");//赋值fM->body_feature_idx
	int a[33] = { 0, 1, 2, 3, 4,5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27, 28, 29, 30, 31, 32 };
	fM->selectbodyfeature.resize(33);
	copy(a, a + 33, fM->selectbodyfeature.begin());
	fM->body_feature_idx1.resize(33);
	fM->body_feature_weights1.resize(33);
	for (int i = 0; i < 33; i++)
	{
		fM->body_feature_idx1[i] = fM->body_feature_idx[a[i]];
		fM->body_feature_weights1[i] = fM->body_feature_weights[a[i]];
	}
	std::vector<int> bodyadd_conset_idx = { 570, 8359, 4152, 5207, 751, 2994, 4093, 677, 4097, 1861 };
	sprintf_s(meshpath, 100, "F:\\scan\\male_without_mouth\\0000.obj");//自己加的
	std::auto_ptr<mns::MyMeshHandler> p_meshHandler(new mns::MyMeshHandler(meshpath, fM->body_feature_idx1, fM->body_feature_weights1));
	Eigen::VectorXd tmpV(fM->nV * 3);
	std::vector<int> smplbody0 = { 3500, 4985, 912, 5333, 1294, 3022, 2967, 6427, 3107, 6530,1619, 1665, 5212, 5165, 4515, 1029, 4486, 1001 };
	std::vector<int> fbmmbody0 = { 4785, 4956, 912, 5306, 1294, 2995, 2940, 6373, 3080, 6476 ,1619, 1665, 5185, 5138, 4488, 1029, 4590, 1001 };
	std::vector<int> smplbody1 = { 3500, 4985, 912, 5333, 1294, 3022, 2967, 6427, 3107, 6530 };
	std::vector<int> fbmmbody1 = { 4785, 4956, 912, 5306, 1294, 2995, 2940, 6373, 3080, 6476 };
	std::vector<int> smplbody2 = { 1619, 1665, 5212, 5165, 4515, 1029, 4486, 1001 };
	std::vector<int> fbmmbody2 = { 1619, 1665, 5185, 5138, 4488, 1029, 4590, 1001 };
	std::vector<int> smplbody3 = { 1970, 1930, 5669, 5560, 6723, 6581, 3323, 3203 };
	std::vector<int> fbmmbody3 = { 1970, 1930, 5642, 5533, 6669, 6527, 3296, 3176 };
	PGMesh * tmpmesh = new PGMesh();
	vector<int> frameset;
	int startframe = 0;
	for (int i = 0; i < (375); i++)
	{
		frameset.push_back(i );
	}
	
	for (int ss = 0; ss < 1; ss++)//for (int ss = 0; ss < frameset.size(); ss++)
	{
		int frame = frameset[ss];
		vector<int> ith_set = { 0 };
		for (int sss = 0; sss < 1; sss++)//for (int sss = 0; sss < ith_set.size(); sss++)
		{
			for (int ii = 1; ii < 101; ii++) {
				int i_th = ith_set[sss];
				int refFrame;
				refFrame = frame;

				//F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\coef\\coef_test_initial_pc.txt
				sprintf_s(meshpath, 100, "E:\\Graphics\\bachelor_thesis\\Code\\smplify-x-master\\output_smpl_server_new\\results\\%02d_img\\coef.txt", ii);//%04d
				readcoefficient(meshpath, x, y);

				//sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\anchorpoint\\%04d.txt", frame);
				//fM->spCM->loadAnchorData(meshpath);
				//fM->spCM->presolve();
				sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\featurepoint\\%04d.txt", frame);
				fM->body_feature_pos1.resize(33);
				readBodyFeaturePos1(meshpath, fM->body_feature_pos1);
				fM->body_feature_pos.resize(33);
				for (int i = 1; i < 33; i++)
				{
					fM->body_feature_pos[i] = fM->body_feature_pos1[i];
				}
				sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\hand_conset\\%04d.txt", frame);
				readHandFeaturePos(meshpath, fM->hand_feature_pos, fM->hand_feature_idx, fM->hand_feature_weights);
				std::vector<std::vector<int>> feet_idx;
				std::vector<std::vector<double>> feet_weights;
				std::vector<Eigen::Vector3d> feet_pos;
				sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\feet_conset\\%04d.txt", frame);
				readHandFeaturePos(meshpath, feet_pos, feet_idx, feet_weights);
				/*sprintf_s(meshpath, 100, "data\\%d\\%d\\faceconset.txt", frame, i_th);
				fM->face_feature_idx.resize(54);
				fM->face_feature_pos.resize(54);
				fM->face_feature_weights.resize(54);
				readFaceFeaturePos(meshpath, fM->face_feature_pos, fM->face_feature_idx, fM->face_feature_weights);*/
				int b[54] = { 9160, 3030, 152, 7173, 9811, 9249, 7258, 8782, 444, 4250, 5242, 1871, 1574, 6445, 4778, 5237, 932, 4336, 4407, 1160, 3098, 1099, 4613, 4578, 3294, 1127, 1050, 870, 4313, 4510, 6551, 3131, 3146, 1216, 2952, 1909, 1976, 6303, 5141, 5401, 5499, 5619, 2794, 1569, 7629, 7690, 4759, 1244, 4167, 6256, 1435, 890, 3063, 4303 };


				std::vector<std::vector<int>> layer1_idx;
				std::vector<std::vector<double>> layer1_weights;
				std::vector<Eigen::Vector3d> layer1_pos;
				sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\layer1_fp\\%04d.txt", frame);
				readHandFeaturePos(meshpath, layer1_pos, layer1_idx, layer1_weights);

				std::vector<std::vector<int>> layer2_idx;//+ feet_idx.size());
				std::vector<std::vector<double>> layer2_weights;// +feet_idx.size());
				std::vector<Eigen::Vector3d> layer2_pos;// +feet_idx.size());
				sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\layer2_fp\\%04d.txt", frame);
				readHandFeaturePos(meshpath, layer2_pos, layer2_idx, layer2_weights);

				std::vector<std::vector<int>> layer3_idx;
				std::vector<std::vector<double>> layer3_weights;
				std::vector<Eigen::Vector3d> layer3_pos;
				sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\layer3_fp\\%04d.txt", frame);
				readHandFeaturePos(meshpath, layer3_pos, layer3_idx, layer3_weights);



				std::vector<int> layer4 = { 15, 16, 17, 18 };
				std::vector<std::vector<int>> layer4_idx(layer4.size());// +facecontour_pos.size());
				std::vector<std::vector<double>> layer4_weights(layer4.size());// +facecontour_pos.size());
				std::vector<Eigen::Vector3d> layer4_pos(layer4.size());

				poseparm = fM->getPoseParm(x, y);
				ofstream outfile1("initialize.txt");
				for (int i = 0; i < 250; i++)
				{
					outfile1 << 0 << endl;
				}
				for (int i = 0; i < 250; i++)
				{
					outfile1 << 0 << endl;
				}
				outfile1.close();
				//cout << poseparm << endl;
				Eigen::MatrixXd V = fM->f_x(x, poseparm);


				Eigen::VectorXd ref = fM->generateShapeWithoutPose(x);
				Eigen::MatrixXd refMatrix(fM->nV, 3);
				for (int i = 0; i < fM->nV; i++)
				{
					refMatrix.row(i) = ref.segment(3 * i, 3);
				}
				fM->spCM->outMesh1(ref, "refmesh.off");
				sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\%04d.ply", frame);//原bodyHands_REGISTRATIONS_B24改为bodyHands_REGISTRATIONS_A08
				OpenMesh::IO::read_mesh(*tmpmesh, meshpath);
				Eigen::MatrixXd smplposMatrix = fM->getVertexMatrix1(tmpmesh);


				Eigen::MatrixXd smplFeaturepos(smplbody0.size(), 3);
				for (int i = 0; i < smplbody0.size(); i++)//加入特征点的位置
				{
					smplFeaturepos.row(i) = smplposMatrix.row(smplbody0[i]);
				}
				Eigen::MatrixXd smplFeaturepos1(smplbody1.size(), 3);
				Eigen::MatrixXd smplFeaturepos2(smplbody2.size(), 3);
				Eigen::MatrixXd smplFeaturepos3(smplbody3.size(), 3);
				std::vector<std::vector<int>> fbmmbody1_idx(smplbody1.size());
				std::vector<std::vector<double>> fbmmbody1_weight(smplbody1.size());
				std::vector<std::vector<int>> fbmmbody2_idx(smplbody2.size());
				std::vector<std::vector<double>> fbmmbody2_weight(smplbody2.size());
				std::vector<std::vector<int>> fbmmbody3_idx(smplbody3.size());
				std::vector<std::vector<double>> fbmmbody3_weight(smplbody3.size());
				for (int i = 0; i < smplbody1.size(); i++)
				{
					smplFeaturepos1.row(i) = smplposMatrix.row(smplbody1[i]);
					std::vector<int> idx = { smplbody1[i] };
					std::vector<double> weight = { 1.0 };
					fbmmbody1_idx[i] = idx;
					fbmmbody1_weight[i] = weight;
				}
				for (int i = 0; i < layer2.size(); i++)
				{
					smplFeaturepos2.row(i) = smplposMatrix.row(smplbody2[i]);
					std::vector<int> idx = { smplbody2[i] };
					std::vector<double> weight = { 1.0 };
					fbmmbody2_idx[i] = idx;
					fbmmbody2_weight[i] = weight;
				}

				for (int i = 0; i < smplbody3.size(); i++)
				{
					smplFeaturepos3.row(i) = smplposMatrix.row(smplbody3[i]);
					std::vector<int> idx = { smplbody1[3] };
					std::vector<double> weight = { 1.0 };
					fbmmbody3_idx[i] = idx;
					fbmmbody3_weight[i] = weight;
				}
				//Eigen::VectorXd currentFeaturepos = fM->getBodyPosfromCurrentV(V, layer1);
				//Eigen::VectorXd currentFeaturepos = fM->getFacePosfromCurrentV(V);
				//优化全局旋转参数R,t
				//fM->findRigidTransform(currentFeaturepos, globalrotate, global_trans, smplFeaturepos);
				//fM->findRigidTransform(currentFeaturepos, globalrotate, global_trans, fM->face_feature_pos);
				Eigen::VectorXd currentFeaturepos = fM->getBodyPosfromCurrentV(V, layer);
				////			//优化全局旋转参数R,t
				fM->findRigidTransform(currentFeaturepos, globalrotate, global_trans, fM->body_feature_pos);
				refMatrix = (globalrotate * refMatrix.transpose()).transpose();
				refMatrix.rowwise() += global_trans.transpose();
				fM->spCM->outMesh1(refMatrix, "refmesh.off");
				fM->spCM->setMesh(refMatrix);
				fM->spCM->calculateFacesNormals();
				fM->spCM->calculateFacesFrame();
				fM->spCM->computeDiEdge();
				Eigen::VectorXd cm0 = fM->spCM->DiEdgeDataMatrix.row(0);

				/*for (int i = 0; i < layer1.size(); i++)
				{
					std::vector<int> idx = fM->body_feature_idx1[layer1[i]];
					std::vector<double> weight = fM->body_feature_weights1[layer1[i]];
					Eigen::Vector3d pos = fM->body_feature_pos[layer1[i]];
					layer1_idx[i] = idx;
					layer1_weights[i] = weight;
					layer1_pos[i] = pos;
				}
				for (int i = 0; i < layer2.size(); i++)
				{
					std::vector<int> idx = fM->body_feature_idx1[layer2[i]];
					std::vector<double> weight = fM->body_feature_weights1[layer2[i]];
					Eigen::Vector3d pos = fM->body_feature_pos[layer2[i]];
					layer2_idx[i] = idx;
					layer2_weights[i] = weight;
					layer2_pos[i] = pos;
				}*/
				/*for (int i = 0; i < feet_idx.size(); i++)
				{
				std::vector<int> idx = feet_idx[i];
				std::vector<double> weight = feet_weights[i];
				Eigen::Vector3d pos = feet_pos[i];
				layer2_idx[layer2.size()+ i] = idx;
				layer2_weights[layer2.size()+i] = weight;
				layer2_pos[layer2.size()+i] = pos;
				}*/
				/*for (int i = 0; i < layer3.size(); i++)
				{
				std::vector<int> idx = fM->body_feature_idx1[layer3[i]];
				std::vector<double> weight = fM->body_feature_weights1[layer3[i]];
				Eigen::Vector3d pos = fM->body_feature_pos[layer3[i]];
				layer3_idx[i] = idx;
				layer3_weights[i] = weight;
				layer3_pos[i] = pos;
				}*/

				for (int i = 0; i < layer4.size(); i++)
				{
					std::vector<int> idx = fM->body_feature_idx1[layer4[i]];
					std::vector<double> weight = fM->body_feature_weights1[layer4[i]];
					Eigen::Vector3d pos = fM->body_feature_pos[layer4[i]];
					layer4_idx[i] = idx;
					layer4_weights[i] = weight;
					layer4_pos[i] = pos;
				}


				Eigen::MatrixXd sRV = globalrotate * (V.transpose());
				sRV.colwise() += global_trans;
				Eigen::MatrixXd sRV1 = sRV.transpose();

				Eigen::MatrixXd sRV2 = sRV1.transpose();
				Eigen::MatrixXd sRV3 = sRV2.transpose();
				sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\init\\%04d.off", frame);
				fM->spCM->outMesh1(sRV3, meshpath);
				sprintf_s(meshpath, 100, "generatedata\\EHF\\%02d_mesh_after_groundtruth.off", ii);
				fM->spCM->outMesh1(sRV1, meshpath);//***输出****refMesh2.off
				p_meshHandler->fitmode = 0;
				//continue; //自己注释掉的
				//PGMesh * pg10 = new PGMesh();//自己注释掉的
				//sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\layer3\\%04d.off", frame-1);//自己注释掉的
				//OpenMesh::IO::read_mesh(*pg10, meshpath);//自己注释掉的
				//Eigen::VectorXd tempV00 = fM->getVertexMatrix(pg10);//自己注释掉的
				//delete pg10;//自己注释掉的
				//Eigen::MatrixXd V10(3, fM->nV);//自己注释掉的
				//for (int i = 0; i < fM->nV; i++)//自己注释掉的
				//{
				//	V10.col(i) = tempV00.segment(3 * i, 3);
				//}


				//after
				//Eigen::MatrixXd V11 = (V10).transpose();//自己注释掉的

				//fM->spCM->outMesh1(V11, "refMesh2.off");//自己注释掉的

				mns::MyMeshHandler* myH = new mns::MyMeshHandler("refMesh2.off", "parameters1_DFAUST.txt", argc, argv);
				//myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\male_expressions\\DE\\C.txt");
				//myH->p_DeformHandler->ReceiveBasis2("C.txt");
				//myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\MFaust\\scripts\\male\\mixedDE\\C.txt");
				myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\female_without_mouth\\DE\\C.txt");
				std::vector<integer> v_id;



				std::vector<int> facefix = { 5298, 1864, 6448, 8387, 4782, 5105, 1651 };
				std::vector<std::vector<int>> facefixidx;
				std::vector<std::vector<double>> facefixweights;
				std::vector<Eigen::Vector3d> facefixpos;
				for (int i = 0; i < facefix.size(); i++)
				{
					vector<int> idx = { facefix[i] };
					vector<double> weight = { 1.0 };
					Eigen::Vector3d coord;
					coord = sRV1.row(facefix[i]).transpose();
					facefixidx.push_back(idx);
					facefixpos.push_back(coord);
					facefixweights.push_back(weight);
				}
				//第一层--脸部
				myH->p_DeformHandler->DeformInit(v_id, V, V);
				myH->p_DeformHandler->setconset(layer1_idx, layer1_weights, layer1_pos);
				myH->p_DeformHandler->addconset(layer2_idx, layer2_weights, layer2_pos);
				//myH->p_DeformHandler->addconset(layer3_idx, layer3_weights, layer3_pos);
				//myH->p_DeformHandler->DeformInit1(v_id, V, V);
				myH->p_DeformHandler->econ_threshold = 1.0;
				myH->p_DeformHandler->Deform(v_id, V, V);
				myH->p_DeformHandler->ModVMat(V);

				sRV1 = V;
				sRV2 = (sRV1.transpose());
				sRV3 = sRV2.transpose();
				sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\layer1\\%04d.off", frame);
				fM->spCM->outMesh1(sRV3, meshpath);

				myH->p_DeformHandler->DeformInit1(v_id, V, V);

				myH->p_DeformHandler->addconset(layer3_idx, layer3_weights, layer3_pos);
				myH->p_DeformHandler->econ_threshold = 1.0;
				myH->p_DeformHandler->Deform(v_id, V, V);
				myH->p_DeformHandler->ModVMat(V);
				sRV1 = V;
				sRV2 = (sRV1.transpose());
				sRV3 = sRV2.transpose();
				sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\layer2\\%04d.off", frame);
				fM->spCM->outMesh1(sRV3, meshpath);

				/* 第二层 身体*/
				PGMesh* pg0 = new PGMesh();
				//outfile1 = ofstream("initialize.txt");
				//for (int i = 0; i < 100; i++)
				//{
				//	outfile1 << 0 << endl;
				//}
				//outfile1.close();
				sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\layer2\\%04d.off", frame);
				OpenMesh::IO::read_mesh(*pg0, meshpath);
				Eigen::VectorXd tempV = fM->getVertexMatrix(pg0);
				delete pg0;
				Eigen::MatrixXd V5(3, fM->nV);
				for (int i = 0; i < fM->nV; i++)
				{
					V5.col(i) = tempV.segment(3 * i, 3);
				}

				/*sRV2 = fromPctosmpl * sRV2;
				sRV = sRV2.transpose();*/
				//after
				Eigen::MatrixXd V6 = (V5).transpose();
				//sRV = sRV2.transpose();

				fM->spCM->outMesh1(V6, "refMesh2.off");
				Eigen::MatrixXd V7(fM->nV, 3);
				myH = new mns::MyMeshHandler("refMesh2.off", "parameters3.txt", argc, argv);
				myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\male_hand\\DE\\newC.txt");
				myH->p_DeformHandler->DeformInit(v_id, V7, V7);
				std::vector<int> handfix = { 4087, 5177, 3038, 2896, 3523, 1678, 1668, 10450, 5221, 4840, 1709 };
				std::vector<std::vector<int>> handfixidx;
				std::vector<std::vector<double>> handfixweights;
				std::vector<Eigen::Vector3d> handfixpos;
				for (int i = 0; i < handfix.size(); i++)
				{
					vector<int> idx = { handfix[i] };
					vector<double> weight = { 1.0 };
					Eigen::Vector3d coord;
					coord = V6.row(handfix[i]).transpose();
					handfixidx.push_back(idx);
					handfixpos.push_back(coord);
					handfixweights.push_back(weight);
				}
				myH->p_DeformHandler->setconset(handfixidx, handfixweights, handfixpos);
				myH->p_DeformHandler->addconset(fM->hand_feature_idx, fM->hand_feature_weights, fM->hand_feature_pos);
				//myH->p_DeformHandler->addconset(layer2_idx, layer2_weights, layer2_pos);
				//myH->p_DeformHandler->DeformInit1(v_id, V, V);
				myH->p_DeformHandler->econ_threshold = 2;
				myH->p_DeformHandler->Deform(v_id, V7, V7);
				myH->p_DeformHandler->ModVMat(V7);

				Eigen::MatrixXd V8 = (V7.transpose());
				Eigen::MatrixXd V9 = V8.transpose();
				sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\layer3\\%04d.off", frame);
				fM->spCM->outMesh1(V9, meshpath);

				//std::vector<int> bodyfix4 = { 3049, 4212, 1460, 4904, 4365, 877, 628, 1394, 5151, 4089, 7591, 1463, 1168, 4496, 4481, 1675, 1314, 4768, 5117 };
				//std::vector<std::vector<int>> bodyfixidx4;
				//std::vector<std::vector<double>> bodyfixweights4;
				//std::vector<Eigen::Vector3d> bodyfixpos4;
				//for (int i = 0; i < bodyfix4.size(); i++)
				//{
				//	vector<int> idx = { bodyfix4[i] };
				//	vector<double> weight = { 1.0 };
				//	Eigen::Vector3d coord;
				//	coord = V7.row(bodyfix4[i]).transpose();
				//	bodyfixidx4.push_back(idx);
				//	bodyfixpos4.push_back(coord);
				//	bodyfixweights4.push_back(weight);
				//}
				//myH->p_DeformHandler->setconset(bodyfixidx, bodyfixweights, bodyfixpos);
				//myH->p_DeformHandler->addconset(layer1_idx, layer1_weights, layer1_pos);
				//myH->p_DeformHandler->addconset(layer2_idx, layer2_weights, layer2_pos);
				//myH->p_DeformHandler->DeformInit1(v_id, V7, V7);
				//myH->p_DeformHandler->Deform(v_id, V7, V7);
				//myH->p_DeformHandler->ModVMat(V7);

				//V8 = s * (fromsmpltoPc * V7.transpose());
				//V8.colwise() += trans;
				//V9 = V8.transpose();
				//sprintf_s(meshpath, 100, "data\\%d\\%d\\layer4.off", frame, i_th);
				//fM->spCM->outMesh1(V9, meshpath);


				////手部姿态
				//PGMesh * pg1 = new PGMesh();
				//outfile1 = ofstream("initialize.txt");
				//for (int i = 0; i < 100; i++)
				//{
				//	outfile1 << 0 << endl;
				//}
				//outfile1.close();
				//sprintf_s(meshpath, 100, "data\\%d\\%d\\layer4.off", frame, i_th);
				//OpenMesh::IO::read_mesh(*pg1, meshpath);
				//Eigen::VectorXd tempV1 = fM->getVertexMatrix(pg1);
				//delete pg1;
				//Eigen::MatrixXd V10(3, fM->nV);
				//for (int i = 0; i < fM->nV; i++)
				//{
				//	V10.col(i) = tempV1.segment(3 * i, 3);
				//}
				//V10.colwise() -= trans;
				//V10 = (1.0 / s) * V10;

				///*sRV2 = fromPctosmpl * sRV2;
				//sRV = sRV2.transpose();*/
				////after
				//Eigen::MatrixXd V11 = (fromPctosmpl * V10).transpose();
				////sRV = sRV2.transpose();

				//fM->spCM->outMesh1(V11, "refMesh2.off");
				//Eigen::MatrixXd V12(fM->nV, 3);
				//myH = new mns::MyMeshHandler("refMesh2.off", "parameters3.txt", argc, argv);
				//myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\male_hand\\DE\\newC.txt");
				//myH->p_DeformHandler->DeformInit(v_id, V12, V12);
				//std::vector<int> handfix = { 4087, 5177, 3038, 2896, 3523, 1678, 1668, 10450, 5221, 4840, 1709 };
				//std::vector<std::vector<int>> handfixidx;
				//std::vector<std::vector<double>> handfixweights;
				//std::vector<Eigen::Vector3d> handfixpos;
				//for (int i = 0; i < handfix.size(); i++)
				//{
				//	vector<int> idx = { handfix[i] };
				//	vector<double> weight = { 1.0 };
				//	Eigen::Vector3d coord;
				//	coord = V11.row(handfix[i]).transpose();
				//	handfixidx.push_back(idx);
				//	handfixpos.push_back(coord);
				//	handfixweights.push_back(weight);
				//}
				//myH->p_DeformHandler->setconset(handfixidx, handfixweights, handfixpos);
				///*myH->p_DeformHandler->setconset(layer1_idx, layer1_weights, layer1_pos);
				//myH->p_DeformHandler->addconset(layer2_idx, layer2_weights, layer2_pos);
				//myH->p_DeformHandler->addconset(layer3_idx, layer3_weights, layer3_pos);*/
				//myH->p_DeformHandler->addconset(fM->hand_feature_idx, fM->hand_feature_weights, fM->hand_feature_pos);

				////myH->p_DeformHandler->DeformInit1(v_id, V, V);
				//myH->p_DeformHandler->Deform(v_id, V12, V12);
				//myH->p_DeformHandler->ModVMat(V12);

				//Eigen::MatrixXd V13 = s * (fromsmpltoPc * V12.transpose());
				//V13.colwise() += trans;
				//Eigen::MatrixXd V14 = V13.transpose();
				//sprintf_s(meshpath, 100, "data\\%d\\%d\\layer5.off", frame, i_th);
				//fM->spCM->outMesh1(V14, meshpath);

			}
		}
	}
	cout << "finish" << endl;
}

void testHumanModel8(int argc, char** argv)
{
	char *meshpath = new char[100];
	sprintf_s(meshpath, 100, "F:\\scan\\male_without_mouth\\0000.obj");
	char *shapepath = new char[100];
	sprintf_s(shapepath, 100, "F:\\scan\\male_without_mouth\\0000.obj");
	char *depath = new char[100];
	sprintf_s(depath, 100, "F:\\scan\\male_without_mouth\\0000.obj");
	PGMesh *pg = new PGMesh();
	OpenMesh::IO::read_mesh(*pg, meshpath);
	Eigen::VectorXd y(69);
	Eigen::VectorXd x(10);
	Eigen::VectorXd facecoef(200);
	Eigen::VectorXd gesturecoef(100);
	Eigen::VectorXd poseparm(500);
	Eigen::Matrix3d globalrotate, fromPctosmpl, fromsmpltoPc;
	double s;
	std::vector<int> layer = { 0, 1, 2, 3, 4, 5,6,7,8,9,10,11,12,13,14,15,16,17,18,19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32 };
	std::vector<int> transvec = { 0, 6, 7, 8, 12, 13, 14 };
	std::vector<int> transvec1 = { 5, 8, 11, 14 };
	std::vector<int> layer1 = { 12, 13, 14, 16, 17, 6, 0, 1, 2, 4, 5, 18, 19, 24, 25, 26, 27, 28, 29, 30, 31, 32 };
	std::vector<int> layer2 = { 20, 21, 22, 23, 7, 8, 10, 11 };
	std::vector<int> layer3 = { 5, 8, 11, 14 };
	Eigen::Vector3d  trans, global_trans;
	//FaceBodyModel * fM = new FaceBodyModel(10, 200, pg);
	FaceBodyModel * fM = new FaceBodyModel(10, 250, pg);
	Eigen::VectorXd out(fM->nV * 3);
	sprintf_s(meshpath, 100, "F:\\scan\\male_expressions\\DE\\newC.txt");
	fM->loadExpressionBases(meshpath, 200);
	/*sprintf_s(meshpath, 100, "F:\\project\\splocs-master\\newC.txt");
	fM->loadPoseBases(meshpath, 400);*/
	sprintf_s(meshpath, 100, "F:\\scan\\MFaust\\scripts\\male\\mixedDE\\newC.txt");
	fM->loadPoseBases(meshpath, 500);
	//sprintf_s(meshpath, 100, "F:\\scan\\female_without_mouth\\DE\\C.txt");
	//fM->loadPoseBases(meshpath, 400);
	sprintf_s(meshpath, 100, "F:\\scan\\male_hand\\DE\\newC.txt");
	fM->loadgestureBases(meshpath, 100);
	sprintf_s(meshpath, 100, "F:\\scan\\smpl\\male_shape\\mean.txt");
	fM->loadMeanShape(meshpath);
	sprintf_s(meshpath, 100, "F:\\scan\\smpl\\male_shape\\shapebasis.txt");
	fM->loadShapeBases(meshpath);

	Eigen::MatrixXd deformedV(fM->nV, 3);
	fM->readBodyFeatureidx1("body_feature_idx1.txt");
	int a[33] = { 0, 1, 2, 3, 4,5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27, 28, 29, 30, 31, 32 };
	fM->selectbodyfeature.resize(33);
	copy(a, a + 33, fM->selectbodyfeature.begin());
	fM->body_feature_idx1.resize(33);
	fM->body_feature_weights1.resize(33);
	for (int i = 0; i < 33; i++)
	{
		fM->body_feature_idx1[i] = fM->body_feature_idx[a[i]];
		fM->body_feature_weights1[i] = fM->body_feature_weights[a[i]];
	}
	std::vector<int> bodyadd_conset_idx = { 570, 8359, 4152, 5207, 751, 2994, 4093, 677, 4097, 1861 };

	std::auto_ptr<mns::MyMeshHandler> p_meshHandler(new mns::MyMeshHandler(meshpath, fM->body_feature_idx1, fM->body_feature_weights1));
	Eigen::VectorXd tmpV(fM->nV * 3);
	std::vector<int> smplbody0 = { 3500, 4985, 912, 5333, 1294, 3022, 2967, 6427, 3107, 6530,1619, 1665, 5212, 5165, 4515, 1029, 4486, 1001 };
	std::vector<int> fbmmbody0 = { 4785, 4956, 912, 5306, 1294, 2995, 2940, 6373, 3080, 6476 ,1619, 1665, 5185, 5138, 4488, 1029, 4590, 1001 };
	std::vector<int> smplbody1 = { 3500, 4985, 912, 5333, 1294, 3022, 2967, 6427, 3107, 6530 };
	std::vector<int> fbmmbody1 = { 4785, 4956, 912, 5306, 1294, 2995, 2940, 6373, 3080, 6476 };
	std::vector<int> smplbody2 = { 1619, 1665, 5212, 5165, 4515, 1029, 4486, 1001 };
	std::vector<int> fbmmbody2 = { 1619, 1665, 5185, 5138, 4488, 1029, 4590, 1001 };
	std::vector<int> smplbody3 = { 1970, 1930, 5669, 5560, 6723, 6581, 3323, 3203 };
	std::vector<int> fbmmbody3 = { 1970, 1930, 5642, 5533, 6669, 6527, 3296, 3176 };
	PGMesh * tmpmesh = new PGMesh();
	vector<int> frameset;
	int startframe = 0;
	for (int i = 200; i < (299); i++)
	{
		frameset.push_back(20 * i);
	}

	for (int ss = 0; ss < frameset.size(); ss++)
	{
		int frame = frameset[ss];
		vector<int> ith_set = { 0 };
		for (int sss = 0; sss < ith_set.size(); sss++)
		{
			int i_th = ith_set[sss];
			int refFrame;
			refFrame = frame;


			sprintf_s(meshpath, 100, "F:\\scan\\MFaust\\scripts\\male\\coef\\%04d.txt", frame);
			readcoefficient(meshpath, x, y);
			sprintf_s(meshpath, 100, "F:\\scan\\MFaust\\scripts\\male\\anchorpoint\\%04d.txt", frame);
			//fM->spCM->loadAnchorData(meshpath);
			//fM->spCM->presolve();
			sprintf_s(meshpath, 100, "F:\\scan\\MFaust\\scripts\\male\\featurepoint\\%04d.txt", frame);
			fM->body_feature_pos1.resize(33);
			readBodyFeaturePos1(meshpath, fM->body_feature_pos1);
			fM->body_feature_pos.resize(33);
			for (int i = 1; i < 33; i++)
			{
				fM->body_feature_pos[i] = fM->body_feature_pos1[i];
			}
			sprintf_s(meshpath, 100, "F:\\scan\\MFaust\\scripts\\male\\hand_conset\\%04d.txt", frame);
			readHandFeaturePos(meshpath, fM->hand_feature_pos, fM->hand_feature_idx, fM->hand_feature_weights);
			std::vector<std::vector<int>> feet_idx;
			std::vector<std::vector<double>> feet_weights;
			std::vector<Eigen::Vector3d> feet_pos;
			sprintf_s(meshpath, 100, "F:\\scan\\MFaust\\scripts\\male\\feet_conset\\%04d.txt", frame);
			readHandFeaturePos(meshpath, feet_pos, feet_idx, feet_weights);
			/*sprintf_s(meshpath, 100, "data\\%d\\%d\\faceconset.txt", frame, i_th);
			fM->face_feature_idx.resize(54);
			fM->face_feature_pos.resize(54);
			fM->face_feature_weights.resize(54);
			readFaceFeaturePos(meshpath, fM->face_feature_pos, fM->face_feature_idx, fM->face_feature_weights);*/
			int b[54] = { 9160, 3030, 152, 7173, 9811, 9249, 7258, 8782, 444, 4250, 5242, 1871, 1574, 6445, 4778, 5237, 932, 4336, 4407, 1160, 3098, 1099, 4613, 4578, 3294, 1127, 1050, 870, 4313, 4510, 6551, 3131, 3146, 1216, 2952, 1909, 1976, 6303, 5141, 5401, 5499, 5619, 2794, 1569, 7629, 7690, 4759, 1244, 4167, 6256, 1435, 890, 3063, 4303 };


			std::vector<std::vector<int>> layer1_idx;
			std::vector<std::vector<double>> layer1_weights;
			std::vector<Eigen::Vector3d> layer1_pos;
			sprintf_s(meshpath, 100, "F:\\scan\\MFaust\\scripts\\male\\layer1_fp\\%04d.txt", frame);
			readHandFeaturePos(meshpath, layer1_pos, layer1_idx, layer1_weights);

			std::vector<std::vector<int>> layer2_idx;//+ feet_idx.size());
			std::vector<std::vector<double>> layer2_weights;// +feet_idx.size());
			std::vector<Eigen::Vector3d> layer2_pos;// +feet_idx.size());
			sprintf_s(meshpath, 100, "F:\\scan\\MFaust\\scripts\\male\\layer2_fp\\%04d.txt", frame);
			readHandFeaturePos(meshpath, layer2_pos, layer2_idx, layer2_weights);

			std::vector<std::vector<int>> layer3_idx;
			std::vector<std::vector<double>> layer3_weights;
			std::vector<Eigen::Vector3d> layer3_pos;
			sprintf_s(meshpath, 100, "F:\\scan\\MFaust\\scripts\\male\\layer3_fp\\%04d.txt", frame);
			readHandFeaturePos(meshpath, layer3_pos, layer3_idx, layer3_weights);



			std::vector<int> layer4 = { 15, 16, 17, 18 };
			std::vector<std::vector<int>> layer4_idx(layer4.size());// +facecontour_pos.size());
			std::vector<std::vector<double>> layer4_weights(layer4.size());// +facecontour_pos.size());
			std::vector<Eigen::Vector3d> layer4_pos(layer4.size());
			poseparm = fM->getPoseParm(x, y);
			ofstream outfile1("initialize.txt");
			for (int i = 0; i < 250; i++)
			{
				outfile1 << 0 << endl;
			}
			for (int i = 0; i < 250; i++)
			{
				outfile1 << 0 << endl;
			}
			outfile1.close();
			//cout << poseparm << endl;
			sprintf_s(shapepath, 100, "F:\\scan\\MFaust\\scripts\\male\\shape\\edition1\\%05d.obj", 2);
			sprintf_s(depath, 100, "F:\\scan\\MFaust\\scripts\\male\\reconst2\\%d.txt", frame);
			Eigen::MatrixXd V = fM->f_x(shapepath, x, poseparm);
			Eigen::VectorXd ref = fM->generateShapeWithoutPose(shapepath);
			Eigen::MatrixXd refMatrix(fM->nV, 3);
			for (int i = 0; i < fM->nV; i++)
			{
				refMatrix.row(i) = ref.segment(3 * i, 3);
			}
			//fM->spCM->outMesh1(ref, "refmesh.off");
			sprintf_s(meshpath, 100, "F:\\scan\\MFaust\\scripts\\male\\%05d.obj", frame);
			OpenMesh::IO::read_mesh(*tmpmesh, meshpath);
			Eigen::MatrixXd smplposMatrix = fM->getVertexMatrix1(tmpmesh);


			Eigen::MatrixXd smplFeaturepos(smplbody0.size(), 3);
			for (int i = 0; i < smplbody0.size(); i++)
			{
				smplFeaturepos.row(i) = smplposMatrix.row(smplbody0[i]);
			}
			Eigen::MatrixXd smplFeaturepos1(smplbody1.size(), 3);
			Eigen::MatrixXd smplFeaturepos2(smplbody2.size(), 3);
			Eigen::MatrixXd smplFeaturepos3(smplbody3.size(), 3);
			std::vector<std::vector<int>> fbmmbody1_idx(smplbody1.size());
			std::vector<std::vector<double>> fbmmbody1_weight(smplbody1.size());
			std::vector<std::vector<int>> fbmmbody2_idx(smplbody2.size());
			std::vector<std::vector<double>> fbmmbody2_weight(smplbody2.size());
			std::vector<std::vector<int>> fbmmbody3_idx(smplbody3.size());
			std::vector<std::vector<double>> fbmmbody3_weight(smplbody3.size());
			for (int i = 0; i < smplbody1.size(); i++)
			{
				smplFeaturepos1.row(i) = smplposMatrix.row(smplbody1[i]);
				std::vector<int> idx = { smplbody1[i] };
				std::vector<double> weight = { 1.0 };
				fbmmbody1_idx[i] = idx;
				fbmmbody1_weight[i] = weight;
			}
			for (int i = 0; i < layer2.size(); i++)
			{
				smplFeaturepos2.row(i) = smplposMatrix.row(smplbody2[i]);
				std::vector<int> idx = { smplbody2[i] };
				std::vector<double> weight = { 1.0 };
				fbmmbody2_idx[i] = idx;
				fbmmbody2_weight[i] = weight;
			}

			for (int i = 0; i < smplbody3.size(); i++)
			{
				smplFeaturepos3.row(i) = smplposMatrix.row(smplbody3[i]);
				std::vector<int> idx = { smplbody1[3] };
				std::vector<double> weight = { 1.0 };
				fbmmbody3_idx[i] = idx;
				fbmmbody3_weight[i] = weight;
			}
			//Eigen::VectorXd currentFeaturepos = fM->getBodyPosfromCurrentV(V, layer1);
			//Eigen::VectorXd currentFeaturepos = fM->getFacePosfromCurrentV(V);
			//优化全局旋转参数R,t
			//fM->findRigidTransform(currentFeaturepos, globalrotate, global_trans, smplFeaturepos);
			//fM->findRigidTransform(currentFeaturepos, globalrotate, global_trans, fM->face_feature_pos);
			Eigen::VectorXd currentFeaturepos = fM->getBodyPosfromCurrentV(V, layer);
			////			//优化全局旋转参数R,t
			fM->findRigidTransform(currentFeaturepos, globalrotate, global_trans, fM->body_feature_pos);
			refMatrix = (globalrotate * refMatrix.transpose()).transpose();
			refMatrix.rowwise() += global_trans.transpose();
			fM->spCM->outMesh1(refMatrix, "refmesh.off");
			fM->spCM->setMesh(refMatrix);
			fM->spCM->calculateFacesNormals();
			fM->spCM->calculateFacesFrame();
			fM->spCM->computeDiEdge();
			Eigen::VectorXd cm0 = fM->spCM->DiEdgeDataMatrix.row(0);

			/*for (int i = 0; i < layer1.size(); i++)
			{
			std::vector<int> idx = fM->body_feature_idx1[layer1[i]];
			std::vector<double> weight = fM->body_feature_weights1[layer1[i]];
			Eigen::Vector3d pos = fM->body_feature_pos[layer1[i]];
			layer1_idx[i] = idx;
			layer1_weights[i] = weight;
			layer1_pos[i] = pos;
			}
			for (int i = 0; i < layer2.size(); i++)
			{
			std::vector<int> idx = fM->body_feature_idx1[layer2[i]];
			std::vector<double> weight = fM->body_feature_weights1[layer2[i]];
			Eigen::Vector3d pos = fM->body_feature_pos[layer2[i]];
			layer2_idx[i] = idx;
			layer2_weights[i] = weight;
			layer2_pos[i] = pos;
			}*/
			/*for (int i = 0; i < feet_idx.size(); i++)
			{
			std::vector<int> idx = feet_idx[i];
			std::vector<double> weight = feet_weights[i];
			Eigen::Vector3d pos = feet_pos[i];
			layer2_idx[layer2.size()+ i] = idx;
			layer2_weights[layer2.size()+i] = weight;
			layer2_pos[layer2.size()+i] = pos;
			}*/
			/*for (int i = 0; i < layer3.size(); i++)
			{
			std::vector<int> idx = fM->body_feature_idx1[layer3[i]];
			std::vector<double> weight = fM->body_feature_weights1[layer3[i]];
			Eigen::Vector3d pos = fM->body_feature_pos[layer3[i]];
			layer3_idx[i] = idx;
			layer3_weights[i] = weight;
			layer3_pos[i] = pos;
			}*/

			for (int i = 0; i < layer4.size(); i++)
			{
				std::vector<int> idx = fM->body_feature_idx1[layer4[i]];
				std::vector<double> weight = fM->body_feature_weights1[layer4[i]];
				Eigen::Vector3d pos = fM->body_feature_pos[layer4[i]];
				layer4_idx[i] = idx;
				layer4_weights[i] = weight;
				layer4_pos[i] = pos;
			}


			Eigen::MatrixXd sRV = globalrotate * (V.transpose());
			sRV.colwise() += global_trans;
			Eigen::MatrixXd sRV1 = sRV.transpose();

			Eigen::MatrixXd sRV2 = sRV1.transpose();
			Eigen::MatrixXd sRV3 = sRV2.transpose();
			sprintf_s(meshpath, 100, "F:\\scan\\MFaust\\scripts\\male\\init_2\\%04d.off", frame);
			fM->spCM->outMesh1(sRV3, meshpath);
			fM->spCM->outMesh1(sRV1, "refMesh2.off");
			p_meshHandler->fitmode = 0;
			continue;
			

			mns::MyMeshHandler* myH = new mns::MyMeshHandler("refMesh2.off", "parameters1_DFAUST.txt", argc, argv);
			//myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\male_expressions\\DE\\C.txt");
			//myH->p_DeformHandler->ReceiveBasis2("C.txt");
			myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\MFaust\\scripts\\male\\mixedDE\\C.txt");
			//myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\female_without_mouth\\DE\\C.txt");
			std::vector<integer> v_id;



			std::vector<int> facefix = { 5298, 1864, 6448, 8387, 4782, 5105, 1651 };
			std::vector<std::vector<int>> facefixidx;
			std::vector<std::vector<double>> facefixweights;
			std::vector<Eigen::Vector3d> facefixpos;
			for (int i = 0; i < facefix.size(); i++)
			{
				vector<int> idx = { facefix[i] };
				vector<double> weight = { 1.0 };
				Eigen::Vector3d coord;
				coord = sRV1.row(facefix[i]).transpose();
				facefixidx.push_back(idx);
				facefixpos.push_back(coord);
				facefixweights.push_back(weight);
			}
			//第一层--脸部
			myH->p_DeformHandler->DeformInit(v_id, V, V);

			myH->p_DeformHandler->setconset(layer1_idx, layer1_weights, layer1_pos);
			myH->p_DeformHandler->addconset(layer2_idx, layer2_weights, layer2_pos);
			//myH->p_DeformHandler->addconset(layer3_idx, layer3_weights, layer3_pos);
			//myH->p_DeformHandler->DeformInit1(v_id, V, V);
			myH->p_DeformHandler->econ_threshold = 1.0;
			myH->p_DeformHandler->Deform(v_id, V, V);
			myH->p_DeformHandler->ModVMat(V);

			sRV1 = V;
			sRV2 = (sRV1.transpose());
			sRV3 = sRV2.transpose();
			sprintf_s(meshpath, 100, "F:\\scan\\MFaust\\scripts\\male\\layer1\\%04d.off", frame);
			fM->spCM->outMesh1(sRV3, meshpath);

			myH->p_DeformHandler->DeformInit1(v_id, V, V);

			myH->p_DeformHandler->addconset(layer3_idx, layer3_weights, layer3_pos);
			myH->p_DeformHandler->econ_threshold = 1.0;
			myH->p_DeformHandler->Deform(v_id, V, V);
			myH->p_DeformHandler->ModVMat(V);
			sRV1 = V;
			sRV2 = (sRV1.transpose());
			sRV3 = sRV2.transpose();
			sprintf_s(meshpath, 100, "F:\\scan\\MFaust\\scripts\\male\\layer2\\%04d.off", frame);
			fM->spCM->outMesh1(sRV3, meshpath);

			/* 第二层 身体*/
			PGMesh * pg0 = new PGMesh();
			outfile1 = ofstream("initialize.txt");
			for (int i = 0; i < 100; i++)
			{
				outfile1 << 0 << endl;
			}
			outfile1.close();
			sprintf_s(meshpath, 100, "F:\\scan\\MFaust\\scripts\\male\\layer2\\%04d.off", frame);
			OpenMesh::IO::read_mesh(*pg0, meshpath);
			Eigen::VectorXd tempV = fM->getVertexMatrix(pg0);
			delete pg0;
			Eigen::MatrixXd V5(3, fM->nV);
			for (int i = 0; i < fM->nV; i++)
			{
				V5.col(i) = tempV.segment(3 * i, 3);
			}

			/*sRV2 = fromPctosmpl * sRV2;
			sRV = sRV2.transpose();*/
			//after
			Eigen::MatrixXd V6 = (V5).transpose();
			//sRV = sRV2.transpose();

			fM->spCM->outMesh1(V6, "refMesh2.off");
			Eigen::MatrixXd V7(fM->nV, 3);
			myH = new mns::MyMeshHandler("refMesh2.off", "parameters3.txt", argc, argv);
			myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\male_hand\\DE\\newC.txt");
			myH->p_DeformHandler->DeformInit(v_id, V7, V7);
			std::vector<int> handfix = { 4087, 5177, 3038, 2896, 3523, 1678, 1668, 10450, 5221, 4840, 1709 };
			std::vector<std::vector<int>> handfixidx;
			std::vector<std::vector<double>> handfixweights;
			std::vector<Eigen::Vector3d> handfixpos;
			for (int i = 0; i < handfix.size(); i++)
			{
				vector<int> idx = { handfix[i] };
				vector<double> weight = { 1.0 };
				Eigen::Vector3d coord;
				coord = V6.row(handfix[i]).transpose();
				handfixidx.push_back(idx);
				handfixpos.push_back(coord);
				handfixweights.push_back(weight);
			}
			myH->p_DeformHandler->setconset(handfixidx, handfixweights, handfixpos);
			myH->p_DeformHandler->addconset(fM->hand_feature_idx, fM->hand_feature_weights, fM->hand_feature_pos);
			//myH->p_DeformHandler->addconset(layer2_idx, layer2_weights, layer2_pos);
			//myH->p_DeformHandler->DeformInit1(v_id, V, V);
			myH->p_DeformHandler->econ_threshold = 2;
			myH->p_DeformHandler->Deform(v_id, V7, V7);
			myH->p_DeformHandler->ModVMat(V7);

			Eigen::MatrixXd V8 = (V7.transpose());
			Eigen::MatrixXd V9 = V8.transpose();
			sprintf_s(meshpath, 100, "F:\\scan\\MFaust\\scripts\\male\\layer3\\%04d.off", frame);
			fM->spCM->outMesh1(V9, meshpath);

			//std::vector<int> bodyfix4 = { 3049, 4212, 1460, 4904, 4365, 877, 628, 1394, 5151, 4089, 7591, 1463, 1168, 4496, 4481, 1675, 1314, 4768, 5117 };
			//std::vector<std::vector<int>> bodyfixidx4;
			//std::vector<std::vector<double>> bodyfixweights4;
			//std::vector<Eigen::Vector3d> bodyfixpos4;
			//for (int i = 0; i < bodyfix4.size(); i++)
			//{
			//	vector<int> idx = { bodyfix4[i] };
			//	vector<double> weight = { 1.0 };
			//	Eigen::Vector3d coord;
			//	coord = V7.row(bodyfix4[i]).transpose();
			//	bodyfixidx4.push_back(idx);
			//	bodyfixpos4.push_back(coord);
			//	bodyfixweights4.push_back(weight);
			//}
			//myH->p_DeformHandler->setconset(bodyfixidx, bodyfixweights, bodyfixpos);
			//myH->p_DeformHandler->addconset(layer1_idx, layer1_weights, layer1_pos);
			//myH->p_DeformHandler->addconset(layer2_idx, layer2_weights, layer2_pos);
			//myH->p_DeformHandler->DeformInit1(v_id, V7, V7);
			//myH->p_DeformHandler->Deform(v_id, V7, V7);
			//myH->p_DeformHandler->ModVMat(V7);

			//V8 = s * (fromsmpltoPc * V7.transpose());
			//V8.colwise() += trans;
			//V9 = V8.transpose();
			//sprintf_s(meshpath, 100, "data\\%d\\%d\\layer4.off", frame, i_th);
			//fM->spCM->outMesh1(V9, meshpath);


			////手部姿态
			//PGMesh * pg1 = new PGMesh();
			//outfile1 = ofstream("initialize.txt");
			//for (int i = 0; i < 100; i++)
			//{
			//	outfile1 << 0 << endl;
			//}
			//outfile1.close();
			//sprintf_s(meshpath, 100, "data\\%d\\%d\\layer4.off", frame, i_th);
			//OpenMesh::IO::read_mesh(*pg1, meshpath);
			//Eigen::VectorXd tempV1 = fM->getVertexMatrix(pg1);
			//delete pg1;
			//Eigen::MatrixXd V10(3, fM->nV);
			//for (int i = 0; i < fM->nV; i++)
			//{
			//	V10.col(i) = tempV1.segment(3 * i, 3);
			//}
			//V10.colwise() -= trans;
			//V10 = (1.0 / s) * V10;

			///*sRV2 = fromPctosmpl * sRV2;
			//sRV = sRV2.transpose();*/
			////after
			//Eigen::MatrixXd V11 = (fromPctosmpl * V10).transpose();
			////sRV = sRV2.transpose();

			//fM->spCM->outMesh1(V11, "refMesh2.off");
			//Eigen::MatrixXd V12(fM->nV, 3);
			//myH = new mns::MyMeshHandler("refMesh2.off", "parameters3.txt", argc, argv);
			//myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\male_hand\\DE\\newC.txt");
			//myH->p_DeformHandler->DeformInit(v_id, V12, V12);
			//std::vector<int> handfix = { 4087, 5177, 3038, 2896, 3523, 1678, 1668, 10450, 5221, 4840, 1709 };
			//std::vector<std::vector<int>> handfixidx;
			//std::vector<std::vector<double>> handfixweights;
			//std::vector<Eigen::Vector3d> handfixpos;
			//for (int i = 0; i < handfix.size(); i++)
			//{
			//	vector<int> idx = { handfix[i] };
			//	vector<double> weight = { 1.0 };
			//	Eigen::Vector3d coord;
			//	coord = V11.row(handfix[i]).transpose();
			//	handfixidx.push_back(idx);
			//	handfixpos.push_back(coord);
			//	handfixweights.push_back(weight);
			//}
			//myH->p_DeformHandler->setconset(handfixidx, handfixweights, handfixpos);
			///*myH->p_DeformHandler->setconset(layer1_idx, layer1_weights, layer1_pos);
			//myH->p_DeformHandler->addconset(layer2_idx, layer2_weights, layer2_pos);
			//myH->p_DeformHandler->addconset(layer3_idx, layer3_weights, layer3_pos);*/
			//myH->p_DeformHandler->addconset(fM->hand_feature_idx, fM->hand_feature_weights, fM->hand_feature_pos);

			////myH->p_DeformHandler->DeformInit1(v_id, V, V);
			//myH->p_DeformHandler->Deform(v_id, V12, V12);
			//myH->p_DeformHandler->ModVMat(V12);

			//Eigen::MatrixXd V13 = s * (fromsmpltoPc * V12.transpose());
			//V13.colwise() += trans;
			//Eigen::MatrixXd V14 = V13.transpose();
			//sprintf_s(meshpath, 100, "data\\%d\\%d\\layer5.off", frame, i_th);
			//fM->spCM->outMesh1(V14, meshpath);


		}
	}

	//cout << y;
}

void testHumanModel8_1(int argc, char** argv)
{
	char *meshpath = new char[100];
	sprintf_s(meshpath, 100, "F:\\scan\\male_without_mouth\\0000.obj");
	char *shapepath = new char[100];
	sprintf_s(shapepath, 100, "F:\\scan\\male_without_mouth\\0000.obj");
	char *depath = new char[100];
	sprintf_s(depath, 100, "F:\\scan\\male_without_mouth\\0000.obj");
	PGMesh *pg = new PGMesh();
	OpenMesh::IO::read_mesh(*pg, meshpath);
	Eigen::VectorXd y(69);
	Eigen::VectorXd x(10);
	Eigen::VectorXd facecoef(200);
	Eigen::VectorXd gesturecoef(100);
	Eigen::VectorXd poseparm(500);
	Eigen::Matrix3d globalrotate, fromPctosmpl, fromsmpltoPc;
	double s;
	std::vector<int> layer = { 0, 1, 2, 3, 4, 5,6,7,8,9,10,11,12,13,14,15,16,17,18,19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32 };
	std::vector<int> transvec = { 0, 6, 7, 8, 12, 13, 14 };
	std::vector<int> transvec1 = { 5, 8, 11, 14 };
	std::vector<int> layer1 = { 12, 13, 14, 16, 17, 6, 0, 1, 2, 4, 5, 18, 19, 24, 25, 26, 27, 28, 29, 30, 31, 32 };
	std::vector<int> layer2 = { 20, 21, 22, 23, 7, 8, 10, 11 };
	std::vector<int> layer3 = { 5, 8, 11, 14 };
	Eigen::Vector3d  trans, global_trans;
	//FaceBodyModel * fM = new FaceBodyModel(10, 200, pg);
	FaceBodyModel * fM = new FaceBodyModel(10, 250, pg);
	Eigen::VectorXd out(fM->nV * 3);
	sprintf_s(meshpath, 100, "F:\\scan\\male_expressions\\DE\\newC.txt");
	fM->loadExpressionBases(meshpath, 200);
	/*sprintf_s(meshpath, 100, "F:\\project\\splocs-master\\newC.txt");
	fM->loadPoseBases(meshpath, 400);*/
	sprintf_s(meshpath, 100, "F:\\scan\\MFaust\\scripts\\male\\mixedDE\\newC.txt");
	fM->loadPoseBases(meshpath, 500);
	//sprintf_s(meshpath, 100, "F:\\scan\\female_without_mouth\\DE\\C.txt");
	//fM->loadPoseBases(meshpath, 400);
	sprintf_s(meshpath, 100, "F:\\scan\\male_hand\\DE\\newC.txt");
	fM->loadgestureBases(meshpath, 100);
	sprintf_s(meshpath, 100, "F:\\scan\\smpl\\male_shape\\mean.txt");
	fM->loadMeanShape(meshpath);
	sprintf_s(meshpath, 100, "F:\\scan\\smpl\\male_shape\\shapebasis.txt");
	fM->loadShapeBases(meshpath);

	Eigen::MatrixXd deformedV(fM->nV, 3);
	fM->readBodyFeatureidx1("body_feature_idx1.txt");
	int a[33] = { 0, 1, 2, 3, 4,5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27, 28, 29, 30, 31, 32 };
	fM->selectbodyfeature.resize(33);
	copy(a, a + 33, fM->selectbodyfeature.begin());
	fM->body_feature_idx1.resize(33);
	fM->body_feature_weights1.resize(33);
	for (int i = 0; i < 33; i++)
	{
		fM->body_feature_idx1[i] = fM->body_feature_idx[a[i]];
		fM->body_feature_weights1[i] = fM->body_feature_weights[a[i]];
	}
	std::vector<int> bodyadd_conset_idx = { 570, 8359, 4152, 5207, 751, 2994, 4093, 677, 4097, 1861 };

	std::auto_ptr<mns::MyMeshHandler> p_meshHandler(new mns::MyMeshHandler(meshpath, fM->body_feature_idx1, fM->body_feature_weights1));
	Eigen::VectorXd tmpV(fM->nV * 3);
	std::vector<int> smplbody0 = { 3500, 4985, 912, 5333, 1294, 3022, 2967, 6427, 3107, 6530,1619, 1665, 5212, 5165, 4515, 1029, 4486, 1001 };
	std::vector<int> fbmmbody0 = { 4785, 4956, 912, 5306, 1294, 2995, 2940, 6373, 3080, 6476 ,1619, 1665, 5185, 5138, 4488, 1029, 4590, 1001 };
	std::vector<int> smplbody1 = { 3500, 4985, 912, 5333, 1294, 3022, 2967, 6427, 3107, 6530 };
	std::vector<int> fbmmbody1 = { 4785, 4956, 912, 5306, 1294, 2995, 2940, 6373, 3080, 6476 };
	std::vector<int> smplbody2 = { 1619, 1665, 5212, 5165, 4515, 1029, 4486, 1001 };
	std::vector<int> fbmmbody2 = { 1619, 1665, 5185, 5138, 4488, 1029, 4590, 1001 };
	std::vector<int> smplbody3 = { 1970, 1930, 5669, 5560, 6723, 6581, 3323, 3203 };
	std::vector<int> fbmmbody3 = { 1970, 1930, 5642, 5533, 6669, 6527, 3296, 3176 };
	PGMesh * tmpmesh = new PGMesh();
	vector<int> frameset;
	int startframe = 0;
	//frameset.push_back(165);
	//frameset.push_back(375);
	//frameset.push_back(37);
	//frameset.push_back(38);
	//frameset.push_back(39);
	
	
	//for (int i = 5; i < (40); i++)
	//{
	//	frameset.push_back(5 * i+1);/*
	//								  frameset.push_back(5 * i + 2);
	//								  frameset.push_back(5 * i + 3);
	//								  frameset.push_back(5 * i + 4);*/
	//}
	//for (int i = 5; i < (40); i++)
	//{
	//	frameset.push_back(5 * i+2);/*
	//							  frameset.push_back(5 * i + 2);
	//							  frameset.push_back(5 * i + 3);
	//							  frameset.push_back(5 * i + 4);*/
	//}
	
	for (int i = 40; i < (80); i++)
	{
		frameset.push_back(5 * i+3);/*
									  frameset.push_back(5 * i + 2);
									  frameset.push_back(5 * i + 3);
									  frameset.push_back(5 * i + 4);*/
	}
	for (int i = 40; i < (80); i++)
	{
		frameset.push_back(5 * i+4);/*
								  frameset.push_back(5 * i + 2);
								  frameset.push_back(5 * i + 3);
								  frameset.push_back(5 * i + 4);*/
	}


	for (int ss = 0; ss < frameset.size(); ss++)
	{
		int frame = frameset[ss];
		vector<int> ith_set = { 0 };
		for (int sss = 0; sss < ith_set.size(); sss++)
		{
			int i_th = ith_set[sss];
			int refFrame;
			refFrame = frame;


			/*sprintf_s(meshpath, 100, "F:\\scan\\MFaust\\scripts\\male\\coef\\%04d.txt", frame);
			readcoefficient(meshpath, x, y);*/
			sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_B24\\anchorpoint\\%04d.txt", frame);
			//fM->spCM->loadAnchorData(meshpath);
			//fM->spCM->presolve();
			sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_B24\\featurepoint\\%04d.txt", frame);
			fM->body_feature_pos1.resize(33);
			readBodyFeaturePos1(meshpath, fM->body_feature_pos1);
			fM->body_feature_pos.resize(33);
			for (int i = 1; i < 33; i++)
			{
				fM->body_feature_pos[i] = fM->body_feature_pos1[i];
			}
			sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_B24\\hand_conset\\%04d.txt", frame);
			readHandFeaturePos(meshpath, fM->hand_feature_pos, fM->hand_feature_idx, fM->hand_feature_weights);
			std::vector<std::vector<int>> feet_idx;
			std::vector<std::vector<double>> feet_weights;
			std::vector<Eigen::Vector3d> feet_pos;
			sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_B24\\feet_conset\\%04d.txt", frame);
			readHandFeaturePos(meshpath, feet_pos, feet_idx, feet_weights);
			/*sprintf_s(meshpath, 100, "data\\%d\\%d\\faceconset.txt", frame, i_th);
			fM->face_feature_idx.resize(54);
			fM->face_feature_pos.resize(54);
			fM->face_feature_weights.resize(54);
			readFaceFeaturePos(meshpath, fM->face_feature_pos, fM->face_feature_idx, fM->face_feature_weights);*/
			int b[54] = { 9160, 3030, 152, 7173, 9811, 9249, 7258, 8782, 444, 4250, 5242, 1871, 1574, 6445, 4778, 5237, 932, 4336, 4407, 1160, 3098, 1099, 4613, 4578, 3294, 1127, 1050, 870, 4313, 4510, 6551, 3131, 3146, 1216, 2952, 1909, 1976, 6303, 5141, 5401, 5499, 5619, 2794, 1569, 7629, 7690, 4759, 1244, 4167, 6256, 1435, 890, 3063, 4303 };


			std::vector<std::vector<int>> layer1_idx;
			std::vector<std::vector<double>> layer1_weights;
			std::vector<Eigen::Vector3d> layer1_pos;
			sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_B24\\layer1_fp\\%04d.txt", frame);
			readHandFeaturePos(meshpath, layer1_pos, layer1_idx, layer1_weights);

			std::vector<std::vector<int>> layer2_idx;//+ feet_idx.size());
			std::vector<std::vector<double>> layer2_weights;// +feet_idx.size());
			std::vector<Eigen::Vector3d> layer2_pos;// +feet_idx.size());
			sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_B24\\layer2_fp\\%04d.txt", frame);
			readHandFeaturePos(meshpath, layer2_pos, layer2_idx, layer2_weights);

			std::vector<std::vector<int>> layer3_idx;
			std::vector<std::vector<double>> layer3_weights;
			std::vector<Eigen::Vector3d> layer3_pos;
			sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_B24\\layer3_fp\\%04d.txt", frame);
			readHandFeaturePos(meshpath, layer3_pos, layer3_idx, layer3_weights);



			std::vector<int> layer4 = { 15, 16, 17, 18 };
			std::vector<std::vector<int>> layer4_idx(layer4.size());// +facecontour_pos.size());
			std::vector<std::vector<double>> layer4_weights(layer4.size());// +facecontour_pos.size());
			std::vector<Eigen::Vector3d> layer4_pos(layer4.size());
			//poseparm = fM->getPoseParm(x, y);
			ofstream outfile1("initialize.txt");
			for (int i = 0; i < 250; i++)
			{
				outfile1 << 0 << endl;
			}
			for (int i = 0; i < 250; i++)
			{
				outfile1 << 0 << endl;
			}
			outfile1.close();
			//cout << poseparm << endl;
			sprintf_s(shapepath, 100, "F:\\scan\\MFaust\\scripts\\male\\shape\\edition1\\%05d.obj", 3);
			sprintf_s(depath, 100, "F:\\scan\\MFaust\\scripts\\male\\reconst2\\%d.txt", frame);
			Eigen::MatrixXd V(fM->nV, 3);
			//Eigen::MatrixXd V = fM->f_x(shapepath, x, poseparm);
			//Eigen::VectorXd ref = fM->generateShapeWithoutPose(shapepath);
			//Eigen::MatrixXd refMatrix(fM->nV, 3);
			//for (int i = 0; i < fM->nV; i++)
			//{
			//	refMatrix.row(i) = ref.segment(3 * i, 3);
			//}
			////fM->spCM->outMesh1(ref, "refmesh.off");
			//sprintf_s(meshpath, 100, "F:\\scan\\MFaust\\scripts\\male\\%05d.obj", frame);
			//OpenMesh::IO::read_mesh(*tmpmesh, meshpath);
			//Eigen::MatrixXd smplposMatrix = fM->getVertexMatrix1(tmpmesh);


			//Eigen::MatrixXd smplFeaturepos(smplbody0.size(), 3);
			//for (int i = 0; i < smplbody0.size(); i++)
			//{
			//	smplFeaturepos.row(i) = smplposMatrix.row(smplbody0[i]);
			//}
			//Eigen::MatrixXd smplFeaturepos1(smplbody1.size(), 3);
			//Eigen::MatrixXd smplFeaturepos2(smplbody2.size(), 3);
			//Eigen::MatrixXd smplFeaturepos3(smplbody3.size(), 3);
			//std::vector<std::vector<int>> fbmmbody1_idx(smplbody1.size());
			//std::vector<std::vector<double>> fbmmbody1_weight(smplbody1.size());
			//std::vector<std::vector<int>> fbmmbody2_idx(smplbody2.size());
			//std::vector<std::vector<double>> fbmmbody2_weight(smplbody2.size());
			//std::vector<std::vector<int>> fbmmbody3_idx(smplbody3.size());
			//std::vector<std::vector<double>> fbmmbody3_weight(smplbody3.size());
			//for (int i = 0; i < smplbody1.size(); i++)
			//{
			//	smplFeaturepos1.row(i) = smplposMatrix.row(smplbody1[i]);
			//	std::vector<int> idx = { smplbody1[i] };
			//	std::vector<double> weight = { 1.0 };
			//	fbmmbody1_idx[i] = idx;
			//	fbmmbody1_weight[i] = weight;
			//}
			//for (int i = 0; i < layer2.size(); i++)
			//{
			//	smplFeaturepos2.row(i) = smplposMatrix.row(smplbody2[i]);
			//	std::vector<int> idx = { smplbody2[i] };
			//	std::vector<double> weight = { 1.0 };
			//	fbmmbody2_idx[i] = idx;
			//	fbmmbody2_weight[i] = weight;
			//}

			//for (int i = 0; i < smplbody3.size(); i++)
			//{
			//	smplFeaturepos3.row(i) = smplposMatrix.row(smplbody3[i]);
			//	std::vector<int> idx = { smplbody1[3] };
			//	std::vector<double> weight = { 1.0 };
			//	fbmmbody3_idx[i] = idx;
			//	fbmmbody3_weight[i] = weight;
			//}
			////Eigen::VectorXd currentFeaturepos = fM->getBodyPosfromCurrentV(V, layer1);
			////Eigen::VectorXd currentFeaturepos = fM->getFacePosfromCurrentV(V);
			////优化全局旋转参数R,t
			////fM->findRigidTransform(currentFeaturepos, globalrotate, global_trans, smplFeaturepos);
			////fM->findRigidTransform(currentFeaturepos, globalrotate, global_trans, fM->face_feature_pos);
			//Eigen::VectorXd currentFeaturepos = fM->getBodyPosfromCurrentV(V, layer);
			//////			//优化全局旋转参数R,t
			//fM->findRigidTransform(currentFeaturepos, globalrotate, global_trans, fM->body_feature_pos);
			//refMatrix = (globalrotate * refMatrix.transpose()).transpose();
			//refMatrix.rowwise() += global_trans.transpose();
			//fM->spCM->outMesh1(refMatrix, "refmesh.off");
			//fM->spCM->setMesh(refMatrix);
			//fM->spCM->calculateFacesNormals();
			//fM->spCM->calculateFacesFrame();
			//fM->spCM->computeDiEdge();
			//Eigen::VectorXd cm0 = fM->spCM->DiEdgeDataMatrix.row(0);

			/*for (int i = 0; i < layer1.size(); i++)
			{
			std::vector<int> idx = fM->body_feature_idx1[layer1[i]];
			std::vector<double> weight = fM->body_feature_weights1[layer1[i]];
			Eigen::Vector3d pos = fM->body_feature_pos[layer1[i]];
			layer1_idx[i] = idx;
			layer1_weights[i] = weight;
			layer1_pos[i] = pos;
			}
			for (int i = 0; i < layer2.size(); i++)
			{
			std::vector<int> idx = fM->body_feature_idx1[layer2[i]];
			std::vector<double> weight = fM->body_feature_weights1[layer2[i]];
			Eigen::Vector3d pos = fM->body_feature_pos[layer2[i]];
			layer2_idx[i] = idx;
			layer2_weights[i] = weight;
			layer2_pos[i] = pos;
			}*/
			/*for (int i = 0; i < feet_idx.size(); i++)
			{
			std::vector<int> idx = feet_idx[i];
			std::vector<double> weight = feet_weights[i];
			Eigen::Vector3d pos = feet_pos[i];
			layer2_idx[layer2.size()+ i] = idx;
			layer2_weights[layer2.size()+i] = weight;
			layer2_pos[layer2.size()+i] = pos;
			}*/
			/*for (int i = 0; i < layer3.size(); i++)
			{
			std::vector<int> idx = fM->body_feature_idx1[layer3[i]];
			std::vector<double> weight = fM->body_feature_weights1[layer3[i]];
			Eigen::Vector3d pos = fM->body_feature_pos[layer3[i]];
			layer3_idx[i] = idx;
			layer3_weights[i] = weight;
			layer3_pos[i] = pos;
			}*/

			/*for (int i = 0; i < layer4.size(); i++)
			{
				std::vector<int> idx = fM->body_feature_idx1[layer4[i]];
				std::vector<double> weight = fM->body_feature_weights1[layer4[i]];
				Eigen::Vector3d pos = fM->body_feature_pos[layer4[i]];
				layer4_idx[i] = idx;
				layer4_weights[i] = weight;
				layer4_pos[i] = pos;
			}


			Eigen::MatrixXd sRV = globalrotate * (V.transpose());
			sRV.colwise() += global_trans;
			Eigen::MatrixXd sRV1 = sRV.transpose();

			Eigen::MatrixXd sRV2 = sRV1.transpose();
			Eigen::MatrixXd sRV3 = sRV2.transpose();
			sprintf_s(meshpath, 100, "F:\\scan\\MFaust\\scripts\\male\\init_2\\%04d.off", frame);
			fM->spCM->outMesh1(sRV3, meshpath);
			fM->spCM->outMesh1(sRV1, "refMesh2.off");
			p_meshHandler->fitmode = 0;
			continue;*/

			sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_B24\\layer2\\%04d.off", frame-1);
			PGMesh * pg11 = new PGMesh();
			outfile1 = ofstream("initialize.txt");
			for (int i = 0; i < 100; i++)
			{
				outfile1 << 0 << endl;
			}
			outfile1.close();
			sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_B24\\layer2\\%04d.off", frame-1);
			OpenMesh::IO::read_mesh(*pg11, meshpath);
			Eigen::VectorXd tempSV1 = fM->getVertexMatrix(pg11);
			delete pg11;
			Eigen::MatrixXd sRV(3, fM->nV);
			for (int i = 0; i < fM->nV; i++)
			{
				sRV.col(i) = tempSV1.segment(3 * i, 3);
			}

			/*sRV2 = fromPctosmpl * sRV2;
			sRV = sRV2.transpose();*/
			//after
			Eigen::MatrixXd sRV1 = (sRV).transpose();
			//sRV = sRV2.transpose();
			fM->spCM->outMesh1(sRV1, meshpath);
			fM->spCM->outMesh1(sRV1, "refMesh2.off");
			mns::MyMeshHandler* myH = new mns::MyMeshHandler("refMesh2.off", "parameters1_DFAUST.txt", argc, argv);
			//myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\male_expressions\\DE\\C.txt");
			//myH->p_DeformHandler->ReceiveBasis2("C.txt");
			myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\MFaust\\scripts\\male\\mixedDE\\C.txt");
			//myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\female_without_mouth\\DE\\C.txt");
			std::vector<integer> v_id;



			std::vector<int> facefix = { 5298, 1864, 6448, 8387, 4782, 5105, 1651 };
			std::vector<std::vector<int>> facefixidx;
			std::vector<std::vector<double>> facefixweights;
			std::vector<Eigen::Vector3d> facefixpos;
			for (int i = 0; i < facefix.size(); i++)
			{
				vector<int> idx = { facefix[i] };
				vector<double> weight = { 1.0 };
				Eigen::Vector3d coord;
				coord = sRV1.row(facefix[i]).transpose();
				facefixidx.push_back(idx);
				facefixpos.push_back(coord);
				facefixweights.push_back(weight);
			}
			//第一层--脸部
			myH->p_DeformHandler->DeformInit(v_id, V, V);

			myH->p_DeformHandler->setconset(layer1_idx, layer1_weights, layer1_pos);
			myH->p_DeformHandler->addconset(layer2_idx, layer2_weights, layer2_pos);
			//myH->p_DeformHandler->addconset(layer3_idx, layer3_weights, layer3_pos);
			//myH->p_DeformHandler->DeformInit1(v_id, V, V);
			myH->p_DeformHandler->econ_threshold = 1.0;
			myH->p_DeformHandler->Deform(v_id, V, V);
			myH->p_DeformHandler->ModVMat(V);

			sRV1 = V;
			sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_B24\\layer1\\%04d.off", frame);
			fM->spCM->outMesh1(sRV1, meshpath);

			myH->p_DeformHandler->DeformInit1(v_id, V, V);

			myH->p_DeformHandler->addconset(layer3_idx, layer3_weights, layer3_pos);
			myH->p_DeformHandler->econ_threshold = 1.0;
			myH->p_DeformHandler->Deform(v_id, V, V);
			myH->p_DeformHandler->ModVMat(V);
			sRV1 = V;
			sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_B24\\layer2\\%04d.off", frame);
			fM->spCM->outMesh1(sRV1, meshpath);

			/* 第二层 身体*/
			PGMesh * pg0 = new PGMesh();
			outfile1 = ofstream("initialize.txt");
			for (int i = 0; i < 100; i++)
			{
				outfile1 << 0 << endl;
			}
			outfile1.close();
			sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_B24\\layer2\\%04d.off", frame);
			OpenMesh::IO::read_mesh(*pg0, meshpath);
			Eigen::VectorXd tempV = fM->getVertexMatrix(pg0);
			delete pg0;
			Eigen::MatrixXd V5(3, fM->nV);
			for (int i = 0; i < fM->nV; i++)
			{
				V5.col(i) = tempV.segment(3 * i, 3);
			}

			/*sRV2 = fromPctosmpl * sRV2;
			sRV = sRV2.transpose();*/
			//after
			Eigen::MatrixXd V6 = (V5).transpose();
			//sRV = sRV2.transpose();

			fM->spCM->outMesh1(V6, "refMesh2.off");
			Eigen::MatrixXd V7(fM->nV, 3);
			myH = new mns::MyMeshHandler("refMesh2.off", "parameters3.txt", argc, argv);
			myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\male_hand\\DE\\newC.txt");
			myH->p_DeformHandler->DeformInit(v_id, V7, V7);
			std::vector<int> handfix = { 4087, 5177, 3038, 2896, 3523, 1678, 1668, 10450, 5221, 4840, 1709 };
			std::vector<std::vector<int>> handfixidx;
			std::vector<std::vector<double>> handfixweights;
			std::vector<Eigen::Vector3d> handfixpos;
			for (int i = 0; i < handfix.size(); i++)
			{
				vector<int> idx = { handfix[i] };
				vector<double> weight = { 1.0 };
				Eigen::Vector3d coord;
				coord = V6.row(handfix[i]).transpose();
				handfixidx.push_back(idx);
				handfixpos.push_back(coord);
				handfixweights.push_back(weight);
			}
			myH->p_DeformHandler->setconset(handfixidx, handfixweights, handfixpos);
			myH->p_DeformHandler->addconset(fM->hand_feature_idx, fM->hand_feature_weights, fM->hand_feature_pos);
			//myH->p_DeformHandler->addconset(layer2_idx, layer2_weights, layer2_pos);
			//myH->p_DeformHandler->DeformInit1(v_id, V, V);
			myH->p_DeformHandler->econ_threshold = 2;
			myH->p_DeformHandler->Deform(v_id, V7, V7);
			myH->p_DeformHandler->ModVMat(V7);

			Eigen::MatrixXd V8 = (V7.transpose());
			Eigen::MatrixXd V9 = V8.transpose();
			sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_B24\\layer3\\%04d.off", frame);
			fM->spCM->outMesh1(V9, meshpath);

			//std::vector<int> bodyfix4 = { 3049, 4212, 1460, 4904, 4365, 877, 628, 1394, 5151, 4089, 7591, 1463, 1168, 4496, 4481, 1675, 1314, 4768, 5117 };
			//std::vector<std::vector<int>> bodyfixidx4;
			//std::vector<std::vector<double>> bodyfixweights4;
			//std::vector<Eigen::Vector3d> bodyfixpos4;
			//for (int i = 0; i < bodyfix4.size(); i++)
			//{
			//	vector<int> idx = { bodyfix4[i] };
			//	vector<double> weight = { 1.0 };
			//	Eigen::Vector3d coord;
			//	coord = V7.row(bodyfix4[i]).transpose();
			//	bodyfixidx4.push_back(idx);
			//	bodyfixpos4.push_back(coord);
			//	bodyfixweights4.push_back(weight);
			//}
			//myH->p_DeformHandler->setconset(bodyfixidx, bodyfixweights, bodyfixpos);
			//myH->p_DeformHandler->addconset(layer1_idx, layer1_weights, layer1_pos);
			//myH->p_DeformHandler->addconset(layer2_idx, layer2_weights, layer2_pos);
			//myH->p_DeformHandler->DeformInit1(v_id, V7, V7);
			//myH->p_DeformHandler->Deform(v_id, V7, V7);
			//myH->p_DeformHandler->ModVMat(V7);

			//V8 = s * (fromsmpltoPc * V7.transpose());
			//V8.colwise() += trans;
			//V9 = V8.transpose();
			//sprintf_s(meshpath, 100, "data\\%d\\%d\\layer4.off", frame, i_th);
			//fM->spCM->outMesh1(V9, meshpath);


			////手部姿态
			//PGMesh * pg1 = new PGMesh();
			//outfile1 = ofstream("initialize.txt");
			//for (int i = 0; i < 100; i++)
			//{
			//	outfile1 << 0 << endl;
			//}
			//outfile1.close();
			//sprintf_s(meshpath, 100, "data\\%d\\%d\\layer4.off", frame, i_th);
			//OpenMesh::IO::read_mesh(*pg1, meshpath);
			//Eigen::VectorXd tempV1 = fM->getVertexMatrix(pg1);
			//delete pg1;
			//Eigen::MatrixXd V10(3, fM->nV);
			//for (int i = 0; i < fM->nV; i++)
			//{
			//	V10.col(i) = tempV1.segment(3 * i, 3);
			//}
			//V10.colwise() -= trans;
			//V10 = (1.0 / s) * V10;

			///*sRV2 = fromPctosmpl * sRV2;
			//sRV = sRV2.transpose();*/
			////after
			//Eigen::MatrixXd V11 = (fromPctosmpl * V10).transpose();
			////sRV = sRV2.transpose();

			//fM->spCM->outMesh1(V11, "refMesh2.off");
			//Eigen::MatrixXd V12(fM->nV, 3);
			//myH = new mns::MyMeshHandler("refMesh2.off", "parameters3.txt", argc, argv);
			//myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\male_hand\\DE\\newC.txt");
			//myH->p_DeformHandler->DeformInit(v_id, V12, V12);
			//std::vector<int> handfix = { 4087, 5177, 3038, 2896, 3523, 1678, 1668, 10450, 5221, 4840, 1709 };
			//std::vector<std::vector<int>> handfixidx;
			//std::vector<std::vector<double>> handfixweights;
			//std::vector<Eigen::Vector3d> handfixpos;
			//for (int i = 0; i < handfix.size(); i++)
			//{
			//	vector<int> idx = { handfix[i] };
			//	vector<double> weight = { 1.0 };
			//	Eigen::Vector3d coord;
			//	coord = V11.row(handfix[i]).transpose();
			//	handfixidx.push_back(idx);
			//	handfixpos.push_back(coord);
			//	handfixweights.push_back(weight);
			//}
			//myH->p_DeformHandler->setconset(handfixidx, handfixweights, handfixpos);
			///*myH->p_DeformHandler->setconset(layer1_idx, layer1_weights, layer1_pos);
			//myH->p_DeformHandler->addconset(layer2_idx, layer2_weights, layer2_pos);
			//myH->p_DeformHandler->addconset(layer3_idx, layer3_weights, layer3_pos);*/
			//myH->p_DeformHandler->addconset(fM->hand_feature_idx, fM->hand_feature_weights, fM->hand_feature_pos);

			////myH->p_DeformHandler->DeformInit1(v_id, V, V);
			//myH->p_DeformHandler->Deform(v_id, V12, V12);
			//myH->p_DeformHandler->ModVMat(V12);

			//Eigen::MatrixXd V13 = s * (fromsmpltoPc * V12.transpose());
			//V13.colwise() += trans;
			//Eigen::MatrixXd V14 = V13.transpose();
			//sprintf_s(meshpath, 100, "data\\%d\\%d\\layer5.off", frame, i_th);
			//fM->spCM->outMesh1(V14, meshpath);


		}
	}

	//cout << y;
}

void testHumanModel8_1_woman(int argc, char** argv)
{
	char *meshpath = new char[100];
	sprintf_s(meshpath, 100, "F:\\scan\\female_without_mouth\\0000.obj");
	char *shapepath = new char[100];
	sprintf_s(shapepath, 100, "F:\\scan\\female_without_mouth\\0000.obj");
	char *depath = new char[100];
	sprintf_s(depath, 100, "F:\\scan\\female_without_mouth\\0000.obj");
	PGMesh *pg = new PGMesh();
	OpenMesh::IO::read_mesh(*pg, meshpath);
	Eigen::VectorXd y(69);
	Eigen::VectorXd x(10);
	Eigen::VectorXd facecoef(200);
	Eigen::VectorXd gesturecoef(100);
	Eigen::VectorXd poseparm(500);
	Eigen::Matrix3d globalrotate, fromPctosmpl, fromsmpltoPc;
	double s;
	std::vector<int> layer = { 0, 1, 2, 3, 4, 5,6,7,8,9,10,11,12,13,14,15,16,17,18,19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32 };
	std::vector<int> transvec = { 0, 6, 7, 8, 12, 13, 14 };
	std::vector<int> transvec1 = { 5, 8, 11, 14 };
	std::vector<int> layer1 = { 12, 13, 14, 16, 17, 6, 0, 1, 2, 4, 5, 18, 19, 24, 25, 26, 27, 28, 29, 30, 31, 32 };
	std::vector<int> layer2 = { 20, 21, 22, 23, 7, 8, 10, 11 };
	std::vector<int> layer3 = { 5, 8, 11, 14 };
	Eigen::Vector3d  trans, global_trans;
	//FaceBodyModel * fM = new FaceBodyModel(10, 200, pg);
	FaceBodyModel * fM = new FaceBodyModel(10, 250, pg);
	Eigen::VectorXd out(fM->nV * 3);
	sprintf_s(meshpath, 100, "F:\\scan\\male_expressions\\DE\\newC.txt");
	fM->loadExpressionBases(meshpath, 200);
	/*sprintf_s(meshpath, 100, "F:\\project\\splocs-master\\newC.txt");
	fM->loadPoseBases(meshpath, 400);*/
	/*sprintf_s(meshpath, 100, "F:\\scan\\MFaust\\scripts\\male\\mixedDE\\newC.txt");
	fM->loadPoseBases(meshpath, 500);*/
	sprintf_s(meshpath, 100, "F:\\scan\\female_without_mouth\\DE\\C.txt");
	fM->loadPoseBases(meshpath, 400);
	sprintf_s(meshpath, 100, "F:\\scan\\male_hand\\DE\\newC.txt");
	fM->loadgestureBases(meshpath, 100);
	sprintf_s(meshpath, 100, "F:\\scan\\smpl\\female_shape\\mean.txt");
	fM->loadMeanShape(meshpath);
	sprintf_s(meshpath, 100, "F:\\scan\\smpl\\female_shape\\shapebasis.txt");
	fM->loadShapeBases(meshpath);

	Eigen::MatrixXd deformedV(fM->nV, 3);
	fM->readBodyFeatureidx1("body_feature_idx1.txt");
	int a[33] = { 0, 1, 2, 3, 4,5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27, 28, 29, 30, 31, 32 };
	fM->selectbodyfeature.resize(33);
	copy(a, a + 33, fM->selectbodyfeature.begin());
	fM->body_feature_idx1.resize(33);
	fM->body_feature_weights1.resize(33);
	for (int i = 0; i < 33; i++)
	{
		fM->body_feature_idx1[i] = fM->body_feature_idx[a[i]];
		fM->body_feature_weights1[i] = fM->body_feature_weights[a[i]];
	}
	std::vector<int> bodyadd_conset_idx = { 570, 8359, 4152, 5207, 751, 2994, 4093, 677, 4097, 1861 };

	std::auto_ptr<mns::MyMeshHandler> p_meshHandler(new mns::MyMeshHandler(meshpath, fM->body_feature_idx1, fM->body_feature_weights1));
	Eigen::VectorXd tmpV(fM->nV * 3);
	std::vector<int> smplbody0 = { 3500, 4985, 912, 5333, 1294, 3022, 2967, 6427, 3107, 6530,1619, 1665, 5212, 5165, 4515, 1029, 4486, 1001 };
	std::vector<int> fbmmbody0 = { 4785, 4956, 912, 5306, 1294, 2995, 2940, 6373, 3080, 6476 ,1619, 1665, 5185, 5138, 4488, 1029, 4590, 1001 };
	std::vector<int> smplbody1 = { 3500, 4985, 912, 5333, 1294, 3022, 2967, 6427, 3107, 6530 };
	std::vector<int> fbmmbody1 = { 4785, 4956, 912, 5306, 1294, 2995, 2940, 6373, 3080, 6476 };
	std::vector<int> smplbody2 = { 1619, 1665, 5212, 5165, 4515, 1029, 4486, 1001 };
	std::vector<int> fbmmbody2 = { 1619, 1665, 5185, 5138, 4488, 1029, 4590, 1001 };
	std::vector<int> smplbody3 = { 1970, 1930, 5669, 5560, 6723, 6581, 3323, 3203 };
	std::vector<int> fbmmbody3 = { 1970, 1930, 5642, 5533, 6669, 6527, 3296, 3176 };
	PGMesh * tmpmesh = new PGMesh();
	vector<int> frameset = {108, 164, 170, 240, 246, 247, 272, 277};
	int startframe = 0;



	/*for (int i = 37; i < (38); i++)
	{
		frameset.push_back(10 * i + 4);
	}*/

	/*for (int i = 0; i < (19); i++)
	{
		frameset.push_back(20 * i + 10);
	}*/

	for (int ss = 0; ss < frameset.size(); ss++)
	{
		int frame = frameset[ss];
		vector<int> ith_set = { 0 };
		for (int sss = 0; sss < ith_set.size(); sss++)
		{
			int i_th = ith_set[sss];
			int refFrame;
			refFrame = frame;


			/*sprintf_s(meshpath, 100, "F:\\scan\\MFaust\\scripts\\male\\coef\\%04d.txt", frame);
			readcoefficient(meshpath, x, y);*/
			sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\anchorpoint\\%04d.txt", frame);
			//fM->spCM->loadAnchorData(meshpath);
			//fM->spCM->presolve();
			sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\featurepoint\\%04d.txt", frame);
			fM->body_feature_pos1.resize(33);
			readBodyFeaturePos1(meshpath, fM->body_feature_pos1);
			fM->body_feature_pos.resize(33);
			for (int i = 1; i < 33; i++)
			{
				fM->body_feature_pos[i] = fM->body_feature_pos1[i];
			}
			sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\hand_conset\\%04d.txt", frame);
			readHandFeaturePos(meshpath, fM->hand_feature_pos, fM->hand_feature_idx, fM->hand_feature_weights);
			std::vector<std::vector<int>> feet_idx;
			std::vector<std::vector<double>> feet_weights;
			std::vector<Eigen::Vector3d> feet_pos;
			sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\feet_conset\\%04d.txt", frame);
			readHandFeaturePos(meshpath, feet_pos, feet_idx, feet_weights);
			/*sprintf_s(meshpath, 100, "data\\%d\\%d\\faceconset.txt", frame, i_th);
			fM->face_feature_idx.resize(54);
			fM->face_feature_pos.resize(54);
			fM->face_feature_weights.resize(54);
			readFaceFeaturePos(meshpath, fM->face_feature_pos, fM->face_feature_idx, fM->face_feature_weights);*/
			int b[54] = { 9160, 3030, 152, 7173, 9811, 9249, 7258, 8782, 444, 4250, 5242, 1871, 1574, 6445, 4778, 5237, 932, 4336, 4407, 1160, 3098, 1099, 4613, 4578, 3294, 1127, 1050, 870, 4313, 4510, 6551, 3131, 3146, 1216, 2952, 1909, 1976, 6303, 5141, 5401, 5499, 5619, 2794, 1569, 7629, 7690, 4759, 1244, 4167, 6256, 1435, 890, 3063, 4303 };


			std::vector<std::vector<int>> layer1_idx;
			std::vector<std::vector<double>> layer1_weights;
			std::vector<Eigen::Vector3d> layer1_pos;
			sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\layer1_fp\\%04d.txt", frame);
			readHandFeaturePos(meshpath, layer1_pos, layer1_idx, layer1_weights);

			std::vector<std::vector<int>> layer2_idx;//+ feet_idx.size());
			std::vector<std::vector<double>> layer2_weights;// +feet_idx.size());
			std::vector<Eigen::Vector3d> layer2_pos;// +feet_idx.size());
			sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\layer2_fp\\%04d.txt", frame);
			readHandFeaturePos(meshpath, layer2_pos, layer2_idx, layer2_weights);

			std::vector<std::vector<int>> layer3_idx;
			std::vector<std::vector<double>> layer3_weights;
			std::vector<Eigen::Vector3d> layer3_pos;
			sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\layer3_fp\\%04d.txt", frame);
			readHandFeaturePos(meshpath, layer3_pos, layer3_idx, layer3_weights);



			std::vector<int> layer4 = { 15, 16, 17, 18 };
			std::vector<std::vector<int>> layer4_idx(layer4.size());// +facecontour_pos.size());
			std::vector<std::vector<double>> layer4_weights(layer4.size());// +facecontour_pos.size());
			std::vector<Eigen::Vector3d> layer4_pos(layer4.size());
			//poseparm = fM->getPoseParm(x, y);
			ofstream outfile1("initialize.txt");
			for (int i = 0; i < 250; i++)
			{
				outfile1 << 0 << endl;
			}
			for (int i = 0; i < 250; i++)
			{
				outfile1 << 0 << endl;
			}
			outfile1.close();
			//cout << poseparm << endl;
			sprintf_s(shapepath, 100, "F:\\scan\\MFaust\\scripts\\male\\shape\\edition1\\%05d.obj", 3);
			sprintf_s(depath, 100, "F:\\scan\\MFaust\\scripts\\male\\reconst2\\%d.txt", frame);
			Eigen::MatrixXd V(fM->nV, 3);
			//Eigen::MatrixXd V = fM->f_x(shapepath, x, poseparm);
			//Eigen::VectorXd ref = fM->generateShapeWithoutPose(shapepath);
			//Eigen::MatrixXd refMatrix(fM->nV, 3);
			//for (int i = 0; i < fM->nV; i++)
			//{
			//	refMatrix.row(i) = ref.segment(3 * i, 3);
			//}
			////fM->spCM->outMesh1(ref, "refmesh.off");
			//sprintf_s(meshpath, 100, "F:\\scan\\MFaust\\scripts\\male\\%05d.obj", frame);
			//OpenMesh::IO::read_mesh(*tmpmesh, meshpath);
			//Eigen::MatrixXd smplposMatrix = fM->getVertexMatrix1(tmpmesh);


			//Eigen::MatrixXd smplFeaturepos(smplbody0.size(), 3);
			//for (int i = 0; i < smplbody0.size(); i++)
			//{
			//	smplFeaturepos.row(i) = smplposMatrix.row(smplbody0[i]);
			//}
			//Eigen::MatrixXd smplFeaturepos1(smplbody1.size(), 3);
			//Eigen::MatrixXd smplFeaturepos2(smplbody2.size(), 3);
			//Eigen::MatrixXd smplFeaturepos3(smplbody3.size(), 3);
			//std::vector<std::vector<int>> fbmmbody1_idx(smplbody1.size());
			//std::vector<std::vector<double>> fbmmbody1_weight(smplbody1.size());
			//std::vector<std::vector<int>> fbmmbody2_idx(smplbody2.size());
			//std::vector<std::vector<double>> fbmmbody2_weight(smplbody2.size());
			//std::vector<std::vector<int>> fbmmbody3_idx(smplbody3.size());
			//std::vector<std::vector<double>> fbmmbody3_weight(smplbody3.size());
			//for (int i = 0; i < smplbody1.size(); i++)
			//{
			//	smplFeaturepos1.row(i) = smplposMatrix.row(smplbody1[i]);
			//	std::vector<int> idx = { smplbody1[i] };
			//	std::vector<double> weight = { 1.0 };
			//	fbmmbody1_idx[i] = idx;
			//	fbmmbody1_weight[i] = weight;
			//}
			//for (int i = 0; i < layer2.size(); i++)
			//{
			//	smplFeaturepos2.row(i) = smplposMatrix.row(smplbody2[i]);
			//	std::vector<int> idx = { smplbody2[i] };
			//	std::vector<double> weight = { 1.0 };
			//	fbmmbody2_idx[i] = idx;
			//	fbmmbody2_weight[i] = weight;
			//}

			//for (int i = 0; i < smplbody3.size(); i++)
			//{
			//	smplFeaturepos3.row(i) = smplposMatrix.row(smplbody3[i]);
			//	std::vector<int> idx = { smplbody1[3] };
			//	std::vector<double> weight = { 1.0 };
			//	fbmmbody3_idx[i] = idx;
			//	fbmmbody3_weight[i] = weight;
			//}
			////Eigen::VectorXd currentFeaturepos = fM->getBodyPosfromCurrentV(V, layer1);
			////Eigen::VectorXd currentFeaturepos = fM->getFacePosfromCurrentV(V);
			////优化全局旋转参数R,t
			////fM->findRigidTransform(currentFeaturepos, globalrotate, global_trans, smplFeaturepos);
			////fM->findRigidTransform(currentFeaturepos, globalrotate, global_trans, fM->face_feature_pos);
			//Eigen::VectorXd currentFeaturepos = fM->getBodyPosfromCurrentV(V, layer);
			//////			//优化全局旋转参数R,t
			//fM->findRigidTransform(currentFeaturepos, globalrotate, global_trans, fM->body_feature_pos);
			//refMatrix = (globalrotate * refMatrix.transpose()).transpose();
			//refMatrix.rowwise() += global_trans.transpose();
			//fM->spCM->outMesh1(refMatrix, "refmesh.off");
			//fM->spCM->setMesh(refMatrix);
			//fM->spCM->calculateFacesNormals();
			//fM->spCM->calculateFacesFrame();
			//fM->spCM->computeDiEdge();
			//Eigen::VectorXd cm0 = fM->spCM->DiEdgeDataMatrix.row(0);

			/*for (int i = 0; i < layer1.size(); i++)
			{
			std::vector<int> idx = fM->body_feature_idx1[layer1[i]];
			std::vector<double> weight = fM->body_feature_weights1[layer1[i]];
			Eigen::Vector3d pos = fM->body_feature_pos[layer1[i]];
			layer1_idx[i] = idx;
			layer1_weights[i] = weight;
			layer1_pos[i] = pos;
			}
			for (int i = 0; i < layer2.size(); i++)
			{
			std::vector<int> idx = fM->body_feature_idx1[layer2[i]];
			std::vector<double> weight = fM->body_feature_weights1[layer2[i]];
			Eigen::Vector3d pos = fM->body_feature_pos[layer2[i]];
			layer2_idx[i] = idx;
			layer2_weights[i] = weight;
			layer2_pos[i] = pos;
			}*/
			/*for (int i = 0; i < feet_idx.size(); i++)
			{
			std::vector<int> idx = feet_idx[i];
			std::vector<double> weight = feet_weights[i];
			Eigen::Vector3d pos = feet_pos[i];
			layer2_idx[layer2.size()+ i] = idx;
			layer2_weights[layer2.size()+i] = weight;
			layer2_pos[layer2.size()+i] = pos;
			}*/
			/*for (int i = 0; i < layer3.size(); i++)
			{
			std::vector<int> idx = fM->body_feature_idx1[layer3[i]];
			std::vector<double> weight = fM->body_feature_weights1[layer3[i]];
			Eigen::Vector3d pos = fM->body_feature_pos[layer3[i]];
			layer3_idx[i] = idx;
			layer3_weights[i] = weight;
			layer3_pos[i] = pos;
			}*/

			/*for (int i = 0; i < layer4.size(); i++)
			{
			std::vector<int> idx = fM->body_feature_idx1[layer4[i]];
			std::vector<double> weight = fM->body_feature_weights1[layer4[i]];
			Eigen::Vector3d pos = fM->body_feature_pos[layer4[i]];
			layer4_idx[i] = idx;
			layer4_weights[i] = weight;
			layer4_pos[i] = pos;
			}


			Eigen::MatrixXd sRV = globalrotate * (V.transpose());
			sRV.colwise() += global_trans;
			Eigen::MatrixXd sRV1 = sRV.transpose();

			Eigen::MatrixXd sRV2 = sRV1.transpose();
			Eigen::MatrixXd sRV3 = sRV2.transpose();
			sprintf_s(meshpath, 100, "F:\\scan\\MFaust\\scripts\\male\\init_2\\%04d.off", frame);
			fM->spCM->outMesh1(sRV3, meshpath);
			fM->spCM->outMesh1(sRV1, "refMesh2.off");
			p_meshHandler->fitmode = 0;
			continue;*/

			sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\layer3\\%04d.off", frame-1);
			PGMesh * pg11 = new PGMesh();
			outfile1 = ofstream("initialize.txt");
			for (int i = 0; i < 100; i++)
			{
				outfile1 << 0 << endl;
			}
			outfile1.close();
			sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\layer3\\%04d.off", frame-1);
			OpenMesh::IO::read_mesh(*pg11, meshpath);
			Eigen::VectorXd tempSV1 = fM->getVertexMatrix(pg11);
			delete pg11;
			Eigen::MatrixXd sRV(3, fM->nV);
			for (int i = 0; i < fM->nV; i++)
			{
				sRV.col(i) = tempSV1.segment(3 * i, 3);
			}

			/*sRV2 = fromPctosmpl * sRV2;
			sRV = sRV2.transpose();*/
			//after
			Eigen::MatrixXd sRV1 = (sRV).transpose();
			//sRV = sRV2.transpose();
			fM->spCM->outMesh1(sRV1, meshpath);
			fM->spCM->outMesh1(sRV1, "refMesh2.off");
			mns::MyMeshHandler* myH = new mns::MyMeshHandler("refMesh2.off", "parameters1_DFAUST.txt", argc, argv);
			//myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\male_expressions\\DE\\C.txt");
			//myH->p_DeformHandler->ReceiveBasis2("C.txt");
			//myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\MFaust\\scripts\\male\\mixedDE\\C.txt");
			myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\female_without_mouth\\DE\\C.txt");
			std::vector<integer> v_id;



			std::vector<int> facefix = { 5298, 1864, 6448, 8387, 4782, 5105, 1651 };
			std::vector<std::vector<int>> facefixidx;
			std::vector<std::vector<double>> facefixweights;
			std::vector<Eigen::Vector3d> facefixpos;
			for (int i = 0; i < facefix.size(); i++)
			{
				vector<int> idx = { facefix[i] };
				vector<double> weight = { 1.0 };
				Eigen::Vector3d coord;
				coord = sRV1.row(facefix[i]).transpose();
				facefixidx.push_back(idx);
				facefixpos.push_back(coord);
				facefixweights.push_back(weight);
			}
			//第一层--脸部
			myH->p_DeformHandler->DeformInit(v_id, V, V);

			myH->p_DeformHandler->setconset(layer1_idx, layer1_weights, layer1_pos);
			myH->p_DeformHandler->addconset(layer2_idx, layer2_weights, layer2_pos);
			//myH->p_DeformHandler->addconset(layer3_idx, layer3_weights, layer3_pos);
			myH->p_DeformHandler->DeformInit1(v_id, V, V);
			myH->p_DeformHandler->econ_threshold = 1.0;
			myH->p_DeformHandler->Deform(v_id, V, V);
			myH->p_DeformHandler->ModVMat(V);

			sRV1 = V;
			sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\layer1\\%04d.off", frame);
			fM->spCM->outMesh1(sRV1, meshpath);

			myH->p_DeformHandler->DeformInit1(v_id, V, V);

			myH->p_DeformHandler->addconset(layer3_idx, layer3_weights, layer3_pos);
			myH->p_DeformHandler->econ_threshold = 1.0;
			myH->p_DeformHandler->Deform(v_id, V, V);
			myH->p_DeformHandler->ModVMat(V);
			sRV1 = V;
			sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\layer2\\%04d.off", frame);
			fM->spCM->outMesh1(sRV1, meshpath);

			/* 第二层 身体*/
			PGMesh * pg0 = new PGMesh();
			outfile1 = ofstream("initialize.txt");
			for (int i = 0; i < 100; i++)
			{
				outfile1 << 0 << endl;
			}
			outfile1.close();
			sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\layer2\\%04d.off", frame);
			OpenMesh::IO::read_mesh(*pg0, meshpath);
			Eigen::VectorXd tempV = fM->getVertexMatrix(pg0);
			delete pg0;
			Eigen::MatrixXd V5(3, fM->nV);
			for (int i = 0; i < fM->nV; i++)
			{
				V5.col(i) = tempV.segment(3 * i, 3);
			}

			/*sRV2 = fromPctosmpl * sRV2;
			sRV = sRV2.transpose();*/
			//after
			Eigen::MatrixXd V6 = (V5).transpose();
			//sRV = sRV2.transpose();

			fM->spCM->outMesh1(V6, "refMesh2.off");
			Eigen::MatrixXd V7(fM->nV, 3);
			myH = new mns::MyMeshHandler("refMesh2.off", "parameters3.txt", argc, argv);
			myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\male_hand\\DE\\newC.txt");
			myH->p_DeformHandler->DeformInit(v_id, V7, V7);
			std::vector<int> handfix = { 4087, 5177, 3038, 2896, 3523, 1678, 1668, 10450, 5221, 4840, 1709 };
			std::vector<std::vector<int>> handfixidx;
			std::vector<std::vector<double>> handfixweights;
			std::vector<Eigen::Vector3d> handfixpos;
			for (int i = 0; i < handfix.size(); i++)
			{
				vector<int> idx = { handfix[i] };
				vector<double> weight = { 1.0 };
				Eigen::Vector3d coord;
				coord = V6.row(handfix[i]).transpose();
				handfixidx.push_back(idx);
				handfixpos.push_back(coord);
				handfixweights.push_back(weight);
			}
			myH->p_DeformHandler->setconset(handfixidx, handfixweights, handfixpos);
			myH->p_DeformHandler->addconset(fM->hand_feature_idx, fM->hand_feature_weights, fM->hand_feature_pos);
			//myH->p_DeformHandler->addconset(layer2_idx, layer2_weights, layer2_pos);
			//myH->p_DeformHandler->DeformInit1(v_id, V, V);
			myH->p_DeformHandler->econ_threshold = 2;
			myH->p_DeformHandler->Deform(v_id, V7, V7);
			myH->p_DeformHandler->ModVMat(V7);

			Eigen::MatrixXd V8 = (V7.transpose());
			Eigen::MatrixXd V9 = V8.transpose();
			sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\layer3\\%04d.off", frame);
			fM->spCM->outMesh1(V9, meshpath);

			//std::vector<int> bodyfix4 = { 3049, 4212, 1460, 4904, 4365, 877, 628, 1394, 5151, 4089, 7591, 1463, 1168, 4496, 4481, 1675, 1314, 4768, 5117 };
			//std::vector<std::vector<int>> bodyfixidx4;
			//std::vector<std::vector<double>> bodyfixweights4;
			//std::vector<Eigen::Vector3d> bodyfixpos4;
			//for (int i = 0; i < bodyfix4.size(); i++)
			//{
			//	vector<int> idx = { bodyfix4[i] };
			//	vector<double> weight = { 1.0 };
			//	Eigen::Vector3d coord;
			//	coord = V7.row(bodyfix4[i]).transpose();
			//	bodyfixidx4.push_back(idx);
			//	bodyfixpos4.push_back(coord);
			//	bodyfixweights4.push_back(weight);
			//}
			//myH->p_DeformHandler->setconset(bodyfixidx, bodyfixweights, bodyfixpos);
			//myH->p_DeformHandler->addconset(layer1_idx, layer1_weights, layer1_pos);
			//myH->p_DeformHandler->addconset(layer2_idx, layer2_weights, layer2_pos);
			//myH->p_DeformHandler->DeformInit1(v_id, V7, V7);
			//myH->p_DeformHandler->Deform(v_id, V7, V7);
			//myH->p_DeformHandler->ModVMat(V7);

			//V8 = s * (fromsmpltoPc * V7.transpose());
			//V8.colwise() += trans;
			//V9 = V8.transpose();
			//sprintf_s(meshpath, 100, "data\\%d\\%d\\layer4.off", frame, i_th);
			//fM->spCM->outMesh1(V9, meshpath);


			////手部姿态
			//PGMesh * pg1 = new PGMesh();
			//outfile1 = ofstream("initialize.txt");
			//for (int i = 0; i < 100; i++)
			//{
			//	outfile1 << 0 << endl;
			//}
			//outfile1.close();
			//sprintf_s(meshpath, 100, "data\\%d\\%d\\layer4.off", frame, i_th);
			//OpenMesh::IO::read_mesh(*pg1, meshpath);
			//Eigen::VectorXd tempV1 = fM->getVertexMatrix(pg1);
			//delete pg1;
			//Eigen::MatrixXd V10(3, fM->nV);
			//for (int i = 0; i < fM->nV; i++)
			//{
			//	V10.col(i) = tempV1.segment(3 * i, 3);
			//}
			//V10.colwise() -= trans;
			//V10 = (1.0 / s) * V10;

			///*sRV2 = fromPctosmpl * sRV2;
			//sRV = sRV2.transpose();*/
			////after
			//Eigen::MatrixXd V11 = (fromPctosmpl * V10).transpose();
			////sRV = sRV2.transpose();

			//fM->spCM->outMesh1(V11, "refMesh2.off");
			//Eigen::MatrixXd V12(fM->nV, 3);
			//myH = new mns::MyMeshHandler("refMesh2.off", "parameters3.txt", argc, argv);
			//myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\male_hand\\DE\\newC.txt");
			//myH->p_DeformHandler->DeformInit(v_id, V12, V12);
			//std::vector<int> handfix = { 4087, 5177, 3038, 2896, 3523, 1678, 1668, 10450, 5221, 4840, 1709 };
			//std::vector<std::vector<int>> handfixidx;
			//std::vector<std::vector<double>> handfixweights;
			//std::vector<Eigen::Vector3d> handfixpos;
			//for (int i = 0; i < handfix.size(); i++)
			//{
			//	vector<int> idx = { handfix[i] };
			//	vector<double> weight = { 1.0 };
			//	Eigen::Vector3d coord;
			//	coord = V11.row(handfix[i]).transpose();
			//	handfixidx.push_back(idx);
			//	handfixpos.push_back(coord);
			//	handfixweights.push_back(weight);
			//}
			//myH->p_DeformHandler->setconset(handfixidx, handfixweights, handfixpos);
			///*myH->p_DeformHandler->setconset(layer1_idx, layer1_weights, layer1_pos);
			//myH->p_DeformHandler->addconset(layer2_idx, layer2_weights, layer2_pos);
			//myH->p_DeformHandler->addconset(layer3_idx, layer3_weights, layer3_pos);*/
			//myH->p_DeformHandler->addconset(fM->hand_feature_idx, fM->hand_feature_weights, fM->hand_feature_pos);

			////myH->p_DeformHandler->DeformInit1(v_id, V, V);
			//myH->p_DeformHandler->Deform(v_id, V12, V12);
			//myH->p_DeformHandler->ModVMat(V12);

			//Eigen::MatrixXd V13 = s * (fromsmpltoPc * V12.transpose());
			//V13.colwise() += trans;
			//Eigen::MatrixXd V14 = V13.transpose();
			//sprintf_s(meshpath, 100, "data\\%d\\%d\\layer5.off", frame, i_th);
			//fM->spCM->outMesh1(V14, meshpath);


		}
	}

	//cout << y;
}

void testHumanModel2(int argc, char** argv)
{
	char *meshpath = new char[100];
	sprintf_s(meshpath, 100, "F:\\scan\\female_without_mouth\\0000.obj");
	PGMesh *pg = new PGMesh();
	OpenMesh::IO::read_mesh(*pg, meshpath);
	Eigen::VectorXd y(69);
	Eigen::VectorXd x(10);
	Eigen::VectorXd facecoef(200);
	Eigen::VectorXd gesturecoef(100);
	Eigen::VectorXd poseparm(400);
	Eigen::Matrix3d globalrotate, fromPctosmpl, fromsmpltoPc;
	double s;
	std::vector<int> layer = { 0, 1, 2, 3, 4, 5,6,7,8,9,10,11,12,13,14,15,16,17,18,19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32 };
	std::vector<int> transvec = { 0, 6, 7, 8, 12, 13, 14 };
	std::vector<int> transvec1 = { 5, 8, 11, 14 };
	std::vector<int> layer1 = { 12, 13, 14, 16, 17, 6, 0, 1, 2, 4, 5, 18, 19, 24, 25, 26, 27, 28, 29, 30, 31, 32 };
	std::vector<int> layer2 = { 20, 21, 22, 23, 7, 8, 10, 11 };
	std::vector<int> layer3 = { 5, 8, 11, 14 };
	Eigen::Vector3d  trans, global_trans;
	FaceBodyModel * fM = new FaceBodyModel(10, 200, pg);
	//FaceBodyModel * fM = new FaceBodyModel(10, 250, pg);
	Eigen::VectorXd out(fM->nV * 3);
	sprintf_s(meshpath, 100, "F:\\scan\\male_expressions\\DE\\newC.txt");
	fM->loadExpressionBases(meshpath, 200);
	/*sprintf_s(meshpath, 100, "F:\\project\\splocs-master\\newC.txt");
	fM->loadPoseBases(meshpath, 400);*/
	/*sprintf_s(meshpath, 100, "F:\\scan\\MFaust\\scripts\\male\\mixedDE\\newC.txt");
	fM->loadPoseBases(meshpath, 500);*/
	sprintf_s(meshpath, 100, "F:\\scan\\female_without_mouth\\DE\\C.txt");
	fM->loadPoseBases(meshpath, 400);
	sprintf_s(meshpath, 100, "F:\\scan\\male_hand\\DE\\newC.txt");
	fM->loadgestureBases(meshpath, 100);
	sprintf_s(meshpath, 100, "F:\\scan\\smpl\\female_shape\\mean.txt");
	fM->loadMeanShape(meshpath);
	sprintf_s(meshpath, 100, "F:\\scan\\smpl\\female_shape\\shapebasis.txt");
	fM->loadShapeBases(meshpath);

	Eigen::MatrixXd deformedV(fM->nV, 3);
	fM->readBodyFeatureidx1("body_feature_idx1.txt");
	int a[33] = {0, 1, 2, 3, 4,5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27, 28, 29, 30, 31, 32 };
	fM->selectbodyfeature.resize(33);
	copy(a, a + 33, fM->selectbodyfeature.begin());
	fM->body_feature_idx1.resize(33);
	fM->body_feature_weights1.resize(33);
	for (int i = 0; i < 33; i++)
	{
		fM->body_feature_idx1[i] = fM->body_feature_idx[a[i]];
		fM->body_feature_weights1[i] = fM->body_feature_weights[a[i]];
	}
	std::vector<int> bodyadd_conset_idx = { 570, 8359, 4152, 5207, 751, 2994, 4093, 677, 4097, 1861 };

	std::auto_ptr<mns::MyMeshHandler> p_meshHandler(new mns::MyMeshHandler(meshpath, fM->body_feature_idx1, fM->body_feature_weights1));
	Eigen::VectorXd tmpV(fM->nV * 3);
	std::vector<int> smplbody0 = { 3500, 4985, 912, 5333, 1294, 3022, 2967, 6427, 3107, 6530,1619, 1665, 5212, 5165, 4515, 1029, 4486, 1001 };
	std::vector<int> fbmmbody0 = { 4785, 4956, 912, 5306, 1294, 2995, 2940, 6373, 3080, 6476 ,1619, 1665, 5185, 5138, 4488, 1029, 4590, 1001 };
	std::vector<int> smplbody1 = { 3500, 4985, 912, 5333, 1294, 3022, 2967, 6427, 3107, 6530 };
	std::vector<int> fbmmbody1 = { 4785, 4956, 912, 5306, 1294, 2995, 2940, 6373, 3080, 6476 };
	std::vector<int> smplbody2 = { 1619, 1665, 5212, 5165, 4515, 1029, 4486, 1001 };
	std::vector<int> fbmmbody2 = { 1619, 1665, 5185, 5138, 4488, 1029, 4590, 1001 };
	std::vector<int> smplbody3 = { 1970, 1930, 5669, 5560, 6723, 6581, 3323, 3203 };
	std::vector<int> fbmmbody3 = { 1970, 1930, 5642, 5533, 6669, 6527, 3296, 3176 };
	PGMesh * tmpmesh = new PGMesh();
	vector<int> frameset ;
	int startframe = 0;
	for (int i = 0; i < (50); i++)
	{
		frameset.push_back( i);
	}
	for (int ss = 0; ss < frameset.size(); ss++)
	{
		int frame = frameset[ss];
		vector<int> ith_set = { 0 };
		for (int sss = 0; sss < ith_set.size(); sss++)
		{
			int i_th = ith_set[sss];
			sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\coef\\%04d.txt", frame);
			readcoefficient(meshpath, x, y);
			sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\anchorpoint\\%04d.txt", frame);
			//fM->spCM->loadAnchorData(meshpath);
			//fM->spCM->presolve();
			sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\featurepoint\\%04d.txt", frame);
			fM->body_feature_pos1.resize(33);
			readBodyFeaturePos1(meshpath, fM->body_feature_pos1);
			fM->body_feature_pos.resize(33);
			for (int i = 1; i < 33; i++)
			{
				fM->body_feature_pos[i] = fM->body_feature_pos1[i];
			}
			sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\hand_conset\\%04d.txt", frame);
			readHandFeaturePos(meshpath, fM->hand_feature_pos, fM->hand_feature_idx, fM->hand_feature_weights);
			std::vector<std::vector<int>> feet_idx;
			std::vector<std::vector<double>> feet_weights;
			std::vector<Eigen::Vector3d> feet_pos;
			sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\feet_conset\\%04d.txt", frame);
			readHandFeaturePos(meshpath, feet_pos, feet_idx, feet_weights);
			/*sprintf_s(meshpath, 100, "data\\%d\\%d\\faceconset.txt", frame, i_th);
			fM->face_feature_idx.resize(54);
			fM->face_feature_pos.resize(54);
			fM->face_feature_weights.resize(54);
			readFaceFeaturePos(meshpath, fM->face_feature_pos, fM->face_feature_idx, fM->face_feature_weights);*/
			int b[54] = { 9160, 3030, 152, 7173, 9811, 9249, 7258, 8782, 444, 4250, 5242, 1871, 1574, 6445, 4778, 5237, 932, 4336, 4407, 1160, 3098, 1099, 4613, 4578, 3294, 1127, 1050, 870, 4313, 4510, 6551, 3131, 3146, 1216, 2952, 1909, 1976, 6303, 5141, 5401, 5499, 5619, 2794, 1569, 7629, 7690, 4759, 1244, 4167, 6256, 1435, 890, 3063, 4303 };
			

			std::vector<std::vector<int>> layer1_idx(layer1.size());
			std::vector<std::vector<double>> layer1_weights(layer1.size());
			std::vector<Eigen::Vector3d> layer1_pos(layer1.size());

			std::vector<std::vector<int>> layer2_idx(layer2.size());//+ feet_idx.size());
			std::vector<std::vector<double>> layer2_weights(layer2.size());// +feet_idx.size());
			std::vector<Eigen::Vector3d> layer2_pos(layer2.size());// +feet_idx.size());

			std::vector<std::vector<int>> layer3_idx = feet_idx;
			std::vector<std::vector<double>> layer3_weights= feet_weights;
			std::vector<Eigen::Vector3d> layer3_pos= feet_pos;



			std::vector<int> layer4 = { 15, 16, 17, 18 };
			std::vector<std::vector<int>> layer4_idx(layer4.size());// +facecontour_pos.size());
			std::vector<std::vector<double>> layer4_weights(layer4.size());// +facecontour_pos.size());
			std::vector<Eigen::Vector3d> layer4_pos(layer4.size());

			poseparm = fM->getPoseParm(x, y);
			ofstream outfile1("initialize.txt");
			for (int i = 0; i < 250; i++)
			{
				outfile1 << 0 << endl;
			}
			for (int i = 0; i < 250; i++)
			{
				outfile1 << 0 << endl;
			}
			outfile1.close();
			//cout << poseparm << endl;
			Eigen::MatrixXd V = fM->f_x(x, poseparm);
			Eigen::VectorXd ref = fM->generateShapeWithoutPose(x);
			Eigen::MatrixXd refMatrix(fM->nV, 3);
			for (int i = 0; i < fM->nV; i++)
			{
				refMatrix.row(i) = ref.segment(3 * i, 3);
			}
			//fM->spCM->outMesh1(ref, "refmesh.off");
			sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\deformed\\%04d.obj", frame);
			OpenMesh::IO::read_mesh(*tmpmesh, meshpath);
			Eigen::MatrixXd smplposMatrix = fM->getVertexMatrix1(tmpmesh);
			

			Eigen::MatrixXd smplFeaturepos(smplbody0.size(), 3);
			for (int i = 0; i < smplbody0.size(); i++)
			{
				smplFeaturepos.row(i) = smplposMatrix.row(smplbody0[i]);
			}
			Eigen::MatrixXd smplFeaturepos1(smplbody1.size(), 3);
			Eigen::MatrixXd smplFeaturepos2(smplbody2.size(), 3);
			Eigen::MatrixXd smplFeaturepos3(smplbody3.size(), 3);
			std::vector<std::vector<int>> fbmmbody1_idx(smplbody1.size());
			std::vector<std::vector<double>> fbmmbody1_weight(smplbody1.size());
			std::vector<std::vector<int>> fbmmbody2_idx(smplbody2.size());
			std::vector<std::vector<double>> fbmmbody2_weight(smplbody2.size());
			std::vector<std::vector<int>> fbmmbody3_idx(smplbody3.size());
			std::vector<std::vector<double>> fbmmbody3_weight(smplbody3.size());
			for (int i = 0; i < smplbody1.size(); i++)
			{
				smplFeaturepos1.row(i) = smplposMatrix.row(smplbody1[i]);
				std::vector<int> idx = { smplbody1[i] };
				std::vector<double> weight = { 1.0 };
				fbmmbody1_idx[i] = idx;
				fbmmbody1_weight[i] = weight;
			}
			for (int i = 0; i < layer2.size(); i++)
			{
				smplFeaturepos2.row(i) = smplposMatrix.row(smplbody2[i]);
				std::vector<int> idx = { smplbody2[i] };
				std::vector<double> weight = { 1.0 };
				fbmmbody2_idx[i] = idx;
				fbmmbody2_weight[i] = weight;
			}

			for (int i = 0; i < smplbody3.size(); i++)
			{
				smplFeaturepos3.row(i) = smplposMatrix.row(smplbody3[i]);
				std::vector<int> idx = { smplbody1[3] };
				std::vector<double> weight = { 1.0 };
				fbmmbody3_idx[i] = idx;
				fbmmbody3_weight[i] = weight;
			}
			//Eigen::VectorXd currentFeaturepos = fM->getBodyPosfromCurrentV(V, layer1);
			//Eigen::VectorXd currentFeaturepos = fM->getFacePosfromCurrentV(V);
			//优化全局旋转参数R,t
			//fM->findRigidTransform(currentFeaturepos, globalrotate, global_trans, smplFeaturepos);
			//fM->findRigidTransform(currentFeaturepos, globalrotate, global_trans, fM->face_feature_pos);
			Eigen::VectorXd currentFeaturepos = fM->getBodyPosfromCurrentV(V, layer);
			////			//优化全局旋转参数R,t
			fM->findRigidTransform(currentFeaturepos, globalrotate, global_trans, fM->body_feature_pos);
			refMatrix = (globalrotate * refMatrix.transpose()).transpose();
			refMatrix.rowwise() += global_trans.transpose();
			fM->spCM->outMesh1(refMatrix, "refmesh.off");
			fM->spCM->setMesh(refMatrix);
			fM->spCM->calculateFacesNormals();
			fM->spCM->calculateFacesFrame();
			fM->spCM->computeDiEdge();
			Eigen::VectorXd cm0 = fM->spCM->DiEdgeDataMatrix.row(0);

			for (int i = 0; i < layer1.size(); i++)
			{
				std::vector<int> idx = fM->body_feature_idx1[layer1[i]];
				std::vector<double> weight = fM->body_feature_weights1[layer1[i]];
				Eigen::Vector3d pos = fM->body_feature_pos[layer1[i]];
				layer1_idx[i] = idx;
				layer1_weights[i] = weight;
				layer1_pos[i] = pos;
			}
			for (int i = 0; i < layer2.size(); i++)
			{
				std::vector<int> idx = fM->body_feature_idx1[layer2[i]];
				std::vector<double> weight = fM->body_feature_weights1[layer2[i]];
				Eigen::Vector3d pos = fM->body_feature_pos[layer2[i]];
				layer2_idx[i] = idx;
				layer2_weights[i] = weight;
				layer2_pos[i] = pos;
			}
			/*for (int i = 0; i < feet_idx.size(); i++)
			{
				std::vector<int> idx = feet_idx[i];
				std::vector<double> weight = feet_weights[i];
				Eigen::Vector3d pos = feet_pos[i];
				layer2_idx[layer2.size()+ i] = idx;
				layer2_weights[layer2.size()+i] = weight;
				layer2_pos[layer2.size()+i] = pos;
			}*/
			/*for (int i = 0; i < layer3.size(); i++)
			{
				std::vector<int> idx = fM->body_feature_idx1[layer3[i]];
				std::vector<double> weight = fM->body_feature_weights1[layer3[i]];
				Eigen::Vector3d pos = fM->body_feature_pos[layer3[i]];
				layer3_idx[i] = idx;
				layer3_weights[i] = weight;
				layer3_pos[i] = pos;
			}*/

			for (int i = 0; i < layer4.size(); i++)
			{
				std::vector<int> idx = fM->body_feature_idx1[layer4[i]];
				std::vector<double> weight = fM->body_feature_weights1[layer4[i]];
				Eigen::Vector3d pos = fM->body_feature_pos[layer4[i]];
				layer4_idx[i] = idx;
				layer4_weights[i] = weight;
				layer4_pos[i] = pos;
			}


			Eigen::MatrixXd sRV = globalrotate * (V.transpose());
			sRV.colwise() += global_trans;
			Eigen::MatrixXd sRV1 = sRV.transpose();

			Eigen::MatrixXd sRV2 = sRV1.transpose();
			Eigen::MatrixXd sRV3 = sRV2.transpose();
			sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\init\\%04d.off", frame);
			fM->spCM->outMesh1(sRV3, meshpath);
			fM->spCM->outMesh1(sRV1, "refMesh2.off");
			p_meshHandler->fitmode = 0;

			mns::MyMeshHandler* myH = new mns::MyMeshHandler("refMesh2.off", "parameters1_DFAUST.txt", argc, argv);
			//myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\male_expressions\\DE\\C.txt");
			//myH->p_DeformHandler->ReceiveBasis2("C.txt");
			//myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\MFaust\\scripts\\male\\mixedDE\\C.txt");
			myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\female_without_mouth\\DE\\C.txt");
			std::vector<integer> v_id;

			

			std::vector<int> facefix = { 5298, 1864, 6448, 8387, 4782, 5105, 1651 };
			std::vector<std::vector<int>> facefixidx;
			std::vector<std::vector<double>> facefixweights;
			std::vector<Eigen::Vector3d> facefixpos;
			for (int i = 0; i < facefix.size(); i++)
			{
				vector<int> idx = { facefix[i] };
				vector<double> weight = { 1.0 };
				Eigen::Vector3d coord;
				coord = sRV1.row(facefix[i]).transpose();
				facefixidx.push_back(idx);
				facefixpos.push_back(coord);
				facefixweights.push_back(weight);
			}
			//第一层--脸部
			myH->p_DeformHandler->DeformInit(v_id, V, V);
			
			myH->p_DeformHandler->setconset(layer1_idx, layer1_weights, layer1_pos);
			myH->p_DeformHandler->addconset(layer2_idx, layer2_weights, layer2_pos);
			//myH->p_DeformHandler->addconset(layer3_idx, layer3_weights, layer3_pos);
			//myH->p_DeformHandler->DeformInit1(v_id, V, V);
			myH->p_DeformHandler->econ_threshold = 1.0;
			myH->p_DeformHandler->Deform(v_id, V, V);
			myH->p_DeformHandler->ModVMat(V);

			sRV1 = V;
			sRV2 = (sRV1.transpose());
			sRV3 = sRV2.transpose();
			sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\layer1\\%04d.off", frame);
			fM->spCM->outMesh1(sRV3, meshpath);

			myH->p_DeformHandler->DeformInit1(v_id, V, V);
			
			myH->p_DeformHandler->addconset(layer3_idx, layer3_weights, layer3_pos);
			myH->p_DeformHandler->econ_threshold = 1.0;
			myH->p_DeformHandler->Deform(v_id, V, V);
			myH->p_DeformHandler->ModVMat(V);
			sRV1 = V;
			sRV2 = (sRV1.transpose());
			sRV3 = sRV2.transpose();
			sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\layer2\\%04d.off", frame);
			fM->spCM->outMesh1(sRV3, meshpath);

			/* 第二层 身体*/
			PGMesh * pg0 = new PGMesh();
			outfile1 = ofstream("initialize.txt");
			for (int i = 0; i < 100; i++)
			{
				outfile1 << 0 << endl;
			}
			outfile1.close();
			sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\layer2\\%04d.off", frame);
			OpenMesh::IO::read_mesh(*pg0, meshpath);
			Eigen::VectorXd tempV = fM->getVertexMatrix(pg0);
			delete pg0;
			Eigen::MatrixXd V5(3, fM->nV);
			for (int i = 0; i < fM->nV; i++)
			{
				V5.col(i) = tempV.segment(3 * i, 3);
			}

			/*sRV2 = fromPctosmpl * sRV2;
			sRV = sRV2.transpose();*/
			//after
			Eigen::MatrixXd V6 = (V5).transpose();
			//sRV = sRV2.transpose();

			fM->spCM->outMesh1(V6, "refMesh2.off");
			Eigen::MatrixXd V7(fM->nV, 3);
			myH = new mns::MyMeshHandler("refMesh2.off", "parameters3.txt", argc, argv);
			myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\male_hand\\DE\\newC.txt");
			myH->p_DeformHandler->DeformInit(v_id, V7, V7);
			std::vector<int> handfix = { 4087, 5177, 3038, 2896, 3523, 1678, 1668, 10450, 5221, 4840, 1709 };
			std::vector<std::vector<int>> handfixidx;
			std::vector<std::vector<double>> handfixweights;
			std::vector<Eigen::Vector3d> handfixpos;
			for (int i = 0; i < handfix.size(); i++)
			{
				vector<int> idx = { handfix[i] };
				vector<double> weight = { 1.0 };
				Eigen::Vector3d coord;
				coord = V6.row(handfix[i]).transpose();
				handfixidx.push_back(idx);
				handfixpos.push_back(coord);
				handfixweights.push_back(weight);
			}
			myH->p_DeformHandler->setconset(handfixidx, handfixweights, handfixpos);
			myH->p_DeformHandler->addconset(fM->hand_feature_idx, fM->hand_feature_weights, fM->hand_feature_pos);
			//myH->p_DeformHandler->addconset(layer2_idx, layer2_weights, layer2_pos);
			//myH->p_DeformHandler->DeformInit1(v_id, V, V);
			myH->p_DeformHandler->econ_threshold = 2;
			myH->p_DeformHandler->Deform(v_id, V7, V7);
			myH->p_DeformHandler->ModVMat(V7);

			Eigen::MatrixXd V8 = ( V7.transpose());
			Eigen::MatrixXd V9 = V8.transpose();
			sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\layer3\\%04d.off", frame);
			fM->spCM->outMesh1(V9, meshpath);

			//std::vector<int> bodyfix4 = { 3049, 4212, 1460, 4904, 4365, 877, 628, 1394, 5151, 4089, 7591, 1463, 1168, 4496, 4481, 1675, 1314, 4768, 5117 };
			//std::vector<std::vector<int>> bodyfixidx4;
			//std::vector<std::vector<double>> bodyfixweights4;
			//std::vector<Eigen::Vector3d> bodyfixpos4;
			//for (int i = 0; i < bodyfix4.size(); i++)
			//{
			//	vector<int> idx = { bodyfix4[i] };
			//	vector<double> weight = { 1.0 };
			//	Eigen::Vector3d coord;
			//	coord = V7.row(bodyfix4[i]).transpose();
			//	bodyfixidx4.push_back(idx);
			//	bodyfixpos4.push_back(coord);
			//	bodyfixweights4.push_back(weight);
			//}
			//myH->p_DeformHandler->setconset(bodyfixidx, bodyfixweights, bodyfixpos);
			//myH->p_DeformHandler->addconset(layer1_idx, layer1_weights, layer1_pos);
			//myH->p_DeformHandler->addconset(layer2_idx, layer2_weights, layer2_pos);
			//myH->p_DeformHandler->DeformInit1(v_id, V7, V7);
			//myH->p_DeformHandler->Deform(v_id, V7, V7);
			//myH->p_DeformHandler->ModVMat(V7);

			//V8 = s * (fromsmpltoPc * V7.transpose());
			//V8.colwise() += trans;
			//V9 = V8.transpose();
			//sprintf_s(meshpath, 100, "data\\%d\\%d\\layer4.off", frame, i_th);
			//fM->spCM->outMesh1(V9, meshpath);


			////手部姿态
			//PGMesh * pg1 = new PGMesh();
			//outfile1 = ofstream("initialize.txt");
			//for (int i = 0; i < 100; i++)
			//{
			//	outfile1 << 0 << endl;
			//}
			//outfile1.close();
			//sprintf_s(meshpath, 100, "data\\%d\\%d\\layer4.off", frame, i_th);
			//OpenMesh::IO::read_mesh(*pg1, meshpath);
			//Eigen::VectorXd tempV1 = fM->getVertexMatrix(pg1);
			//delete pg1;
			//Eigen::MatrixXd V10(3, fM->nV);
			//for (int i = 0; i < fM->nV; i++)
			//{
			//	V10.col(i) = tempV1.segment(3 * i, 3);
			//}
			//V10.colwise() -= trans;
			//V10 = (1.0 / s) * V10;

			///*sRV2 = fromPctosmpl * sRV2;
			//sRV = sRV2.transpose();*/
			////after
			//Eigen::MatrixXd V11 = (fromPctosmpl * V10).transpose();
			////sRV = sRV2.transpose();

			//fM->spCM->outMesh1(V11, "refMesh2.off");
			//Eigen::MatrixXd V12(fM->nV, 3);
			//myH = new mns::MyMeshHandler("refMesh2.off", "parameters3.txt", argc, argv);
			//myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\male_hand\\DE\\newC.txt");
			//myH->p_DeformHandler->DeformInit(v_id, V12, V12);
			//std::vector<int> handfix = { 4087, 5177, 3038, 2896, 3523, 1678, 1668, 10450, 5221, 4840, 1709 };
			//std::vector<std::vector<int>> handfixidx;
			//std::vector<std::vector<double>> handfixweights;
			//std::vector<Eigen::Vector3d> handfixpos;
			//for (int i = 0; i < handfix.size(); i++)
			//{
			//	vector<int> idx = { handfix[i] };
			//	vector<double> weight = { 1.0 };
			//	Eigen::Vector3d coord;
			//	coord = V11.row(handfix[i]).transpose();
			//	handfixidx.push_back(idx);
			//	handfixpos.push_back(coord);
			//	handfixweights.push_back(weight);
			//}
			//myH->p_DeformHandler->setconset(handfixidx, handfixweights, handfixpos);
			///*myH->p_DeformHandler->setconset(layer1_idx, layer1_weights, layer1_pos);
			//myH->p_DeformHandler->addconset(layer2_idx, layer2_weights, layer2_pos);
			//myH->p_DeformHandler->addconset(layer3_idx, layer3_weights, layer3_pos);*/
			//myH->p_DeformHandler->addconset(fM->hand_feature_idx, fM->hand_feature_weights, fM->hand_feature_pos);

			////myH->p_DeformHandler->DeformInit1(v_id, V, V);
			//myH->p_DeformHandler->Deform(v_id, V12, V12);
			//myH->p_DeformHandler->ModVMat(V12);

			//Eigen::MatrixXd V13 = s * (fromsmpltoPc * V12.transpose());
			//V13.colwise() += trans;
			//Eigen::MatrixXd V14 = V13.transpose();
			//sprintf_s(meshpath, 100, "data\\%d\\%d\\layer5.off", frame, i_th);
			//fM->spCM->outMesh1(V14, meshpath);


		}
	}

	//cout << y;
}


void testHumanModel1(int argc, char** argv)
{
	clock_t start, finish;
	char *meshpath = new char[100];
	sprintf_s(meshpath, 100, "F:\\scan\\male_without_mouth\\0000.obj");
	PGMesh *pg = new PGMesh();
	OpenMesh::IO::read_mesh(*pg, meshpath);
	Eigen::VectorXd y(69);
	Eigen::VectorXd x(10);
	Eigen::VectorXd facecoef(200);
	Eigen::VectorXd gesturecoef(100);
	Eigen::VectorXd poseparm(400);
	Eigen::Matrix3d globalrotate, fromPctosmpl, fromsmpltoPc;
	double s;
	std::vector<int> layer = { 0, 1, 2, 3, 4, 5,6,7,8,9,10,11,12,13,14,15,16,17,18,19, 20, 21, 22, 23, 24, 25, 26, 27, 28 };
	std::vector<int> transvec = { 0, 6, 7, 8, 12, 13, 14};
	std::vector<int> transvec1 = {5, 8, 11, 14 };
	std::vector<int> layer1 = {  3, 6, 9, 12 };
	std::vector<int> layer2 = { 4, 7, 10, 13, 5, 8, 11, 14 };
	std::vector<int> layer3 = { 5, 8, 11, 14 };
	Eigen::Vector3d  trans, global_trans;
	FaceBodyModel * fM = new FaceBodyModel(10, 200, pg);
	Eigen::VectorXd out(fM->nV * 3);
	sprintf_s(meshpath, 100, "F:\\scan\\male_expressions\\DE\\newC.txt");
	fM->loadExpressionBases(meshpath, 200);
	//sprintf_s(meshpath, 100, "F:\\project\\splocs-master\\newC.txt");
	sprintf_s(meshpath, 100, "F:\\scan\\MFaust\\scripts\\male\\mixedDE\\newC.txt");
	//sprintf_s(meshpath, 100, "F:\\scan\\female_without_mouth\\DE\\C.txt");
	fM->loadPoseBases(meshpath, 500);
	sprintf_s(meshpath, 100, "F:\\scan\\male_hand\\DE\\newC.txt");
	fM->loadgestureBases(meshpath, 100);
	sprintf_s(meshpath, 100, "F:\\scan\\smpl\\male_shape\\mean.txt");
	fM->loadMeanShape(meshpath);
	sprintf_s(meshpath, 100, "F:\\scan\\smpl\\male_shape\\shapebasis.txt");
	fM->loadShapeBases(meshpath);

	//output shape basis
	//for (int i = 0; i < 10; i++)
	//{
	//	Eigen::VectorXd s_coef(10);
	//	s_coef.setZero();
	//	s_coef(i) = 3.0;
	//	Eigen::VectorXd y_out = fM->generateShapeWithoutPose(s_coef);
	//	std::vector<double> error(fM->nV);
	//	for (int j = 0; j < fM->nV; j++)
	//	{
	//		error[j] = (fM->C_S.block<1, 3>(i, 3 * j + 0)) * (fM->C_S.block<1, 3>(i, 3 * j + 0).transpose());
	//		error[j] = sqrt(error[j]);
	//	}
	//	sprintf_s(meshpath, 100, "shape_component\\mesh_%04d.obj", i);
	//	fM->spCM->colorize(error, y_out, meshpath);
	//}

	Eigen::MatrixXd deformedV(fM->nV, 3);
	fM->readBodyFeatureidx("body_feature_idx.txt");
	int a[19] = { 12, 15, 0, 16, 18, 20, 1, 4, 7, 17, 19, 21, 2, 5, 8, 24, 25, 26, 27 };
	fM->selectbodyfeature.resize(19);
	copy(a, a + 19, fM->selectbodyfeature.begin());
	fM->body_feature_idx1.resize(19);
	fM->body_feature_weights1.resize(19);
	for (int i = 0; i < 19; i++)
	{
		fM->body_feature_idx1[i] = fM->body_feature_idx[a[i]];
		fM->body_feature_weights1[i] = fM->body_feature_weights[a[i]];
	}
	std::vector<int> bodyadd_conset_idx = {570, 8359, 4152, 5207, 751, 2994, 4093, 677, 4097, 1861};
	for (int i = 19; i < 29; i++)
	{
		std::vector<int> idx = { bodyadd_conset_idx [i-19]};
		std::vector<double> weight = {1.0};
	}
	std::auto_ptr<mns::MyMeshHandler> p_meshHandler(new mns::MyMeshHandler(meshpath, fM->body_feature_idx1, fM->body_feature_weights1));
	Eigen::VectorXd tmpV(fM->nV * 3);
	std::vector<int> smplbody0 = { 3500, 4985, 912, 5333, 1294, 3022, 2967, 6427, 3107, 6530,1619, 1665, 5212, 5165, 4515, 1029, 4486, 1001 };
	std::vector<int> fbmmbody0 = { 4785, 4956, 912, 5306, 1294, 2995, 2940, 6373, 3080, 6476 ,1619, 1665, 5185, 5138, 4488, 1029, 4590, 1001 };
	std::vector<int> smplbody1 = {3500, 4985, 912, 5333, 1294, 3022, 2967, 6427, 3107, 6530};
	std::vector<int> fbmmbody1 = {4785, 4956, 912, 5306, 1294, 2995, 2940, 6373, 3080, 6476};
	std::vector<int> smplbody2 = {1619, 1665, 5212, 5165, 4515, 1029, 4486, 1001};
	std::vector<int> fbmmbody2 = {1619, 1665, 5185, 5138, 4488, 1029, 4590, 1001};
	std::vector<int> smplbody3 = {1970, 1930, 5669, 5560, 6723, 6581, 3323, 3203};
	std::vector<int> fbmmbody3 = {1970, 1930, 5642, 5533, 6669, 6527, 3296, 3176};
	PGMesh * tmpmesh = new PGMesh();
	vector<int> frameset ;
	for (int i = 5198; i < 5250; i++)
	{
		frameset.push_back(i);
	}
	ofstream timeRecord("time_statistic.txt");
	for (int ss = 0; ss < frameset.size(); ss++)
	{
		int frame = frameset[ss];
		vector<int> ith_set = {0};
		for (int sss = 0; sss < ith_set.size(); sss++)
		{
			int i_th = ith_set[sss];
			sprintf_s(meshpath, 100, "data\\%d\\%d\\coef.txt", frame, i_th);
			readcoefficient(meshpath, x, y);
			sprintf_s(meshpath, 100, "data\\%d\\%d\\global.txt", frame, i_th);
			read_global(meshpath, globalrotate, fromPctosmpl, fromsmpltoPc, s, trans, global_trans);
			sprintf_s(meshpath, 100, "data\\%d\\%d\\anchorpoint.txt", frame, i_th);
			//fM->spCM->loadAnchorData(meshpath);
			//fM->spCM->presolve();
			sprintf_s(meshpath, 100, "data\\%d\\%d\\bodyconset.txt", frame, i_th);
			fM->body_feature_pos1.resize(19);
			readBodyFeaturePos(meshpath, fM->body_feature_pos1);
			fM->body_feature_pos.resize(19);
			for (int i = 1; i < 19; i++)
			{
				fM->body_feature_pos[i] = fM->body_feature_pos1[i];
			}
			sprintf_s(meshpath, 100, "data\\%d\\%d\\handconset.txt", frame, i_th);
			fM->hand_feature_idx.resize(42);
			fM->hand_feature_weights.resize(42);
			fM->hand_feature_pos.resize(42);
			readHandFeaturePos(meshpath, fM->hand_feature_pos, fM->hand_feature_idx, fM->hand_feature_weights);
			/*for (int h_idx = 0; h_idx < fM->hand_feature_pos.size(); h_idx++)
			{
				Eigen::Vector3d v0 = fM->hand_feature_pos[h_idx];
				Eigen::Vector3d v1 = (1 / s) * (fromPctosmpl * (v0 - trans));
				fM->hand_feature_pos[h_idx] = v1;
			}*/
			sprintf_s(meshpath, 100, "data\\%d\\%d\\faceconset.txt", frame, i_th);
			fM->face_feature_idx.resize(54);
			fM->face_feature_pos.resize(54);
			fM->face_feature_weights.resize(54);
			readHandFeaturePos(meshpath, fM->face_feature_pos, fM->face_feature_idx, fM->face_feature_weights);
			
			int b[54] = { 9160, 3030, 152, 7173, 9811, 9249, 7258, 8782, 444, 4250, 5242, 1871, 1574, 6445, 4778, 5237, 932, 4336, 4407, 1160, 3098, 1099, 4613, 4578, 3294, 1127, 1050, 870, 4313, 4510, 6551, 3131, 3146, 1216, 2952, 1909, 1976, 6303, 5141, 5401, 5499, 5619, 2794, 1569, 7629, 7690, 4759, 1244, 4167, 6256, 1435, 890, 3063, 4303 };
			/*for (int i = 0; i < 54; i++)
			{
			fM->face_feature_idx[54 + i].resize(1);
			fM->face_feature_idx[54 + i][0] = b[i];
			fM->face_feature_idx[54 + i][0] = b[i];
			fM->face_feature_weights[54 + i].resize(1);
			fM->face_feature_weights[54 + i][0] = 1.0;
			}*/
			
			/*std::vector<std::vector<int>> layer1_idx(fbmmbody1.size());
			std::vector<std::vector<double>> layer1_weights(fbmmbody1.size());
			std::vector<Eigen::Vector3d> layer1_pos(fbmmbody1.size());

			std::vector<std::vector<int>> layer2_idx(fbmmbody2.size() + 4);
			std::vector<std::vector<double>> layer2_weights(fbmmbody2.size()+ 4);
			std::vector<Eigen::Vector3d> layer2_pos(fbmmbody2.size()+ 4);
			
			std::vector<std::vector<int>> layer3_idx(fbmmbody3.size());
			std::vector<std::vector<double>> layer3_weights(fbmmbody3.size());
			std::vector<Eigen::Vector3d> layer3_pos(fbmmbody3.size());*/

			std::vector<std::vector<int>> layer1_idx(layer1.size());
			std::vector<std::vector<double>> layer1_weights(layer1.size());
			std::vector<Eigen::Vector3d> layer1_pos(layer1.size());

			std::vector<std::vector<int>> layer2_idx(layer2.size());
			std::vector<std::vector<double>> layer2_weights(layer2.size());
			std::vector<Eigen::Vector3d> layer2_pos(layer2.size() );

			std::vector<std::vector<int>> layer3_idx(layer3.size());
			std::vector<std::vector<double>> layer3_weights(layer3.size());
			std::vector<Eigen::Vector3d> layer3_pos(layer3.size());

			std::vector<std::vector<int>> facecontour_idx;
			std::vector<std::vector<double>> facecontour_weights;
			std::vector<Eigen::Vector3d> facecontour_pos;
			//sprintf_s(meshpath, 100, "data\\%d\\%d\\facecontourconset.txt", frame, i_th);
			//readFacecontourFeaturePos(meshpath, facecontour_pos, facecontour_idx, facecontour_weights);

			std::vector<int> layer4 = { 15, 16, 17, 18 };
			std::vector<std::vector<int>> layer4_idx(layer4.size());// +facecontour_pos.size());
			std::vector<std::vector<double>> layer4_weights(layer4.size());// +facecontour_pos.size());
			std::vector<Eigen::Vector3d> layer4_pos(layer4.size());// +facecontour_pos.size());

			/*for (int i = 0; i < facecontour_idx.size(); i++)
			{
				fM->face_feature_idx.push_back(facecontour_idx[i]);
				fM->face_feature_weights.push_back(facecontour_weights[i]);
				fM->face_feature_pos.push_back(facecontour_pos[i]);
			}*/
			start = clock();
			poseparm = fM->getPoseParm(x, y);
			finish = clock();
			double totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
			timeRecord << totaltime << " ";
			ofstream outfile1("initialize.txt");
			for (int i = 0; i < 200; i++)
			{
				outfile1 << 0 << endl;
			}
			for (int i = 0; i < 200; i++)
			{
				outfile1 << 0 << endl;
			}
			outfile1.close();
			//cout << poseparm << endl;
			start = clock();
			Eigen::MatrixXd V = fM->f_x(x, poseparm);
			finish = clock();
			totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
			timeRecord << totaltime << " ";
			Eigen::VectorXd ref = fM->generateShapeWithoutPose(x);
			Eigen::MatrixXd refMatrix(fM->nV, 3);
			for (int i = 0; i < fM->nV; i++)
			{
				refMatrix.row(i) = ref.segment(3 * i, 3);
			}
			//fM->spCM->outMesh1(ref, "refmesh.off");
			sprintf_s(meshpath, 100, "data\\%d\\%d\\deformed.obj", frame, i_th);
			OpenMesh::IO::read_mesh(*tmpmesh, meshpath);
			Eigen::MatrixXd smplposMatrix = fM->getVertexMatrix1(tmpmesh);
			smplposMatrix.rowwise() -= trans.transpose();
			Eigen::MatrixXd smplposMatrix_T = (1.0 / s) * (smplposMatrix.transpose());
			smplposMatrix_T = fromPctosmpl * (smplposMatrix_T);
			smplposMatrix = (smplposMatrix_T.transpose());
			
			Eigen::MatrixXd smplFeaturepos(smplbody0.size(), 3);
			for (int i = 0; i < smplbody0.size(); i++)
			{
				smplFeaturepos.row(i) = smplposMatrix.row(smplbody0[i]);
			}
			Eigen::MatrixXd smplFeaturepos1(smplbody1.size(), 3);
			Eigen::MatrixXd smplFeaturepos2(smplbody2.size(), 3);
			Eigen::MatrixXd smplFeaturepos3(smplbody3.size(), 3);
			std::vector<std::vector<int>> fbmmbody1_idx(smplbody1.size());
			std::vector<std::vector<double>> fbmmbody1_weight(smplbody1.size());
			std::vector<std::vector<int>> fbmmbody2_idx(smplbody2.size());
			std::vector<std::vector<double>> fbmmbody2_weight(smplbody2.size());
			std::vector<std::vector<int>> fbmmbody3_idx(smplbody3.size());
			std::vector<std::vector<double>> fbmmbody3_weight(smplbody3.size());
			for (int i = 0; i < smplbody1.size(); i++)
			{
				smplFeaturepos1.row(i) = smplposMatrix.row(smplbody1[i]);
				std::vector<int> idx = { smplbody1[i] };
				std::vector<double> weight = { 1.0 };
				fbmmbody1_idx[i] = idx;
				fbmmbody1_weight[i] = weight;
			}
			for (int i = 0; i < smplbody2.size(); i++)
			{
				smplFeaturepos2.row(i) = smplposMatrix.row(smplbody2[i]);
				std::vector<int> idx = { smplbody2[i] };
				std::vector<double> weight = { 1.0 };
				fbmmbody2_idx[i] = idx;
				fbmmbody2_weight[i] = weight;
			}
			for (int i = 0; i < smplbody3.size(); i++)
			{
				smplFeaturepos3.row(i) = smplposMatrix.row(smplbody3[i]);
				std::vector<int> idx = { smplbody1[3] };
				std::vector<double> weight = { 1.0 };
				fbmmbody3_idx[i] = idx;
				fbmmbody3_weight[i] = weight;
			}
			//Eigen::VectorXd currentFeaturepos = fM->getBodyPosfromCurrentV1(V, fbmmbody0);
			Eigen::VectorXd currentFeaturepos = fM->getFacePosfromCurrentV(V);
			//优化全局旋转参数R,t
			//fM->findRigidTransform(currentFeaturepos, globalrotate, global_trans, smplFeaturepos);
			fM->findRigidTransform(currentFeaturepos, globalrotate, global_trans, fM->face_feature_pos);
			//Eigen::VectorXd currentFeaturepos = fM->getBodyPosfromCurrentV(V, layer);
			////			//优化全局旋转参数R,t
			//fM->findRigidTransform(currentFeaturepos, globalrotate, global_trans, fM->body_feature_pos);
			refMatrix = (globalrotate * refMatrix.transpose()).transpose();
			refMatrix.rowwise() += global_trans.transpose();
			fM->spCM->outMesh1(refMatrix, "refmesh.off");
			sprintf_s(meshpath, 100, "data\\%d\\%d\\ref.off", frame, i_th);
			fM->spCM->outMesh1(refMatrix, meshpath);
			fM->spCM->setMesh(refMatrix);
			fM->spCM->calculateFacesNormals();
			fM->spCM->calculateFacesFrame();
			fM->spCM->computeDiEdge();
			Eigen::VectorXd cm0 = fM->spCM->DiEdgeDataMatrix.row(0);
			
			for (int i = 0; i < layer1.size(); i++)
			{
				std::vector<int> idx = fM->body_feature_idx1[layer1[i]];
				std::vector<double> weight = fM->body_feature_weights1[layer1[i]];
				Eigen::Vector3d pos = fM->body_feature_pos[layer1[i]];
				layer1_idx[i] = idx;
				layer1_weights[i] = weight;
				layer1_pos[i] = pos;
			}
			for (int i = 0; i < layer2.size(); i++)
			{
				std::vector<int> idx = fM->body_feature_idx1[layer2[i]];
				std::vector<double> weight = fM->body_feature_weights1[layer2[i]];
				Eigen::Vector3d pos = fM->body_feature_pos[layer2[i]];
				layer2_idx[i] = idx;
				layer2_weights[i] = weight;
				layer2_pos[i] = pos;
			}
			for (int i = 0; i < layer3.size(); i++)
			{
				std::vector<int> idx = fM->body_feature_idx1[layer3[i]];
				std::vector<double> weight = fM->body_feature_weights1[layer3[i]];
				Eigen::Vector3d pos = fM->body_feature_pos[layer3[i]];
				layer3_idx[i] = idx;
				layer3_weights[i] = weight;
				layer3_pos[i] = pos;
			}

			for (int i = 0; i < layer4.size(); i++)
			{
				std::vector<int> idx = fM->body_feature_idx1[layer4[i]];
				std::vector<double> weight = fM->body_feature_weights1[layer4[i]];
				Eigen::Vector3d pos = fM->body_feature_pos[layer4[i]];
				layer4_idx[i] = idx;
				layer4_weights[i] = weight;
				layer4_pos[i] = pos;
			}

			/*for (int i = 0; i < facecontour_idx.size(); i++)
			{
				std::vector<int> idx = (facecontour_idx[i]);
				std::vector<double> weight = (facecontour_weights[i]);
				Eigen::Vector3d pos = (facecontour_pos[i]);
				layer4_idx[layer4.size() + i] = idx;
				layer4_weights[layer4.size()+i] = weight;
				layer4_pos[layer4.size()+ i] = pos;
			}*/


			////
			//Eigen::MatrixXd sRV = globalrotate * (V.transpose());
			//输入表面特征点
			/*for (int i = 0; i < fbmmbody1.size(); i++)
			{
				std::vector<int> idx = {fbmmbody1[i]};
				std::vector<double> weight = { 1.0 };
				Eigen::Vector3d pos = smplposMatrix.row(smplbody1[i]).transpose();
				layer1_idx[i] = idx;
				layer1_weights[i] = weight;
				layer1_pos[i] = pos;
			}
			for (int i = 0; i < fbmmbody2.size(); i++)
			{
				std::vector<int> idx = { fbmmbody2[i] };
				std::vector<double> weight = { 1.0 };
				Eigen::Vector3d pos = smplposMatrix.row(smplbody2[i]).transpose();
				layer2_idx[i] = idx;
				layer2_weights[i] = weight;
				layer2_pos[i] = pos;
			}
			layer2_idx[layer2_idx.size() - 4] = fM->body_feature_idx1[15];
			layer2_idx[layer2_idx.size() - 3] = fM->body_feature_idx1[16];
			layer2_idx[layer2_idx.size() - 2] = fM->body_feature_idx1[17];
			layer2_idx[layer2_idx.size() - 1] = fM->body_feature_idx1[18];
			layer2_weights[layer2_idx.size() - 4] = fM->body_feature_weights1[15];
			layer2_weights[layer2_idx.size() - 3] = fM->body_feature_weights1[16];
			layer2_weights[layer2_idx.size() - 2] = fM->body_feature_weights1[17];
			layer2_weights[layer2_idx.size() - 1] = fM->body_feature_weights1[18];
			layer2_pos[layer2_idx.size() - 4] = fM->body_feature_pos1[15];
			layer2_pos[layer2_idx.size() - 3] = fM->body_feature_pos1[16];
			layer2_pos[layer2_idx.size() - 2] = fM->body_feature_pos1[17];
			layer2_pos[layer2_idx.size() - 1] = fM->body_feature_pos1[18];
			for (int i = 0; i < fbmmbody3.size(); i++)
			{
				std::vector<int> idx = { fbmmbody3[i] };
				std::vector<double> weight = { 1.0 };
				Eigen::Vector3d pos = smplposMatrix.row(smplbody3[i]).transpose();
				layer3_idx[i] = idx;
				layer3_weights[i] = weight;
				layer3_pos[i] = pos;
			}*/
			Eigen::MatrixXd sRV = globalrotate * (V.transpose());
			sRV.colwise() += global_trans;
			Eigen::MatrixXd sRV1 = sRV.transpose();
			//currentFeaturepos = fM->getBodyPosfromCurrentV(sRV1, transvec);

			Eigen::MatrixXd sRV2 = s * (fromsmpltoPc * sRV1.transpose());
			sRV2.colwise() += trans;
			Eigen::MatrixXd sRV3 = sRV2.transpose();
			sprintf_s(meshpath, 100, "data\\%d\\%d\\initial_global.off", frame, i_th);
			fM->spCM->outMesh1(sRV3, meshpath);
			sprintf_s(meshpath, 100, "data\\%d\\%d\\initial_reconstruction.off", frame, i_th);
			fM->spCM->outMesh1(sRV1, meshpath);
			fM->spCM->outMesh1(sRV1, "refMesh2.off"); 
			p_meshHandler->fitmode = 0;

			mns::MyMeshHandler* myH = new mns::MyMeshHandler("refMesh2.off", "parameters2.txt", argc, argv);
			myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\male_expressions\\DE\\C.txt");
			//myH->p_DeformHandler->ReceiveBasis2("C.txt");
			std::vector<integer> v_id;

			//第一层--脸部
			//myH->p_DeformHandler->DeformInit(v_id, V, V);
			//
			//myH->p_DeformHandler->setconset(layer1_idx, layer1_weights, layer1_pos);
			//myH->p_DeformHandler->addconset(layer2_idx, layer2_weights, layer2_pos);
			////myH->p_DeformHandler->DeformInit1(v_id, V, V);
			//myH->p_DeformHandler->econ_threshold = 1.0;
			//myH->p_DeformHandler->Deform(v_id, V, V);
			//myH->p_DeformHandler->ModVMat(V);

			//sRV1 = V;
			//sRV2 = s * (fromsmpltoPc * sRV1.transpose());
			//sRV2.colwise() += trans;
			//sRV3 = sRV2.transpose();
			//sprintf_s(meshpath, 100, "data\\%d\\%d\\layer1.off", frame, i_th);
			//fM->spCM->outMesh1(sRV3, meshpath);
			//fM->spCM->outMesh1(sRV1, "refMesh2.off");


			//std::vector<std::vector<int>> faceidx1;
			//std::vector<std::vector<double>> faceweights1;
			//std::vector<Eigen::Vector3d> facepos1;
			//for (int i = 0; i < 20; i++)
			//{
			//	faceidx1.push_back(fM->face_feature_idx[i]);
			//	faceweights1.push_back(fM->face_feature_weights[i]);
			//	facepos1.push_back(fM->face_feature_pos[i]);
			//}
			//std::vector<std::vector<int>> faceidx2;
			//std::vector<std::vector<double>> faceweights2;
			//std::vector<Eigen::Vector3d> facepos2;
			//for (int i = 20; i < 54; i++)
			//{
			//	faceidx2.push_back(fM->face_feature_idx[i]);
			//	faceweights2.push_back(fM->face_feature_weights[i]);
			//	facepos2.push_back(fM->face_feature_pos[i]);
			//}
			//myH = new mns::MyMeshHandler("refMesh2.off", "parameters2.txt", argc, argv);
			//myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\male_expressions\\DE\\C.txt");
			//myH->p_DeformHandler->DeformInit(v_id, V, V);
			//myH->p_DeformHandler->setconset(facefixidx, facefixweights, facefixpos);
			////myH->p_DeformHandler->addconset(faceidx2, faceweights2, facepos2);
			//myH->p_DeformHandler->addconset(faceidx1, faceweights1, facepos1);
			//myH->p_DeformHandler->econ_threshold = 1.0;
			////myH->p_DeformHandler->Deform(v_id, V, V);
			//myH->p_DeformHandler->ModVMat(V);
			//sRV1 = V;
			//sRV2 = s * (fromsmpltoPc * sRV1.transpose());
			//sRV2.colwise() += trans;
			//sRV3 = sRV2.transpose();
			//sprintf_s(meshpath, 100, "data\\%d\\%d\\layer2.off", frame, i_th);
			//fM->spCM->outMesh1(sRV3, meshpath);


			std::vector<int> facefix = { 5298, 1864, 6448, 8387, 4782, 5105, 1651 };
			std::vector<std::vector<int>> facefixidx;
			std::vector<std::vector<double>> facefixweights;
			std::vector<Eigen::Vector3d> facefixpos;
			for (int i = 0; i < facefix.size(); i++)
			{
				vector<int> idx = { facefix[i] };
				vector<double> weight = { 1.0 };
				Eigen::Vector3d coord;
				coord = sRV1.row(facefix[i]).transpose();
				facefixidx.push_back(idx);
				facefixpos.push_back(coord);
				facefixweights.push_back(weight);
			}
			fM->face_feature_idx[23][0] = 3673;
			fM->face_feature_idx[24][0] = 8926;
			fM->face_feature_idx[25][0] = 9087;
			fM->face_feature_idx[26][0] = 9585;
			fM->face_feature_idx[27][0] = 3719;
			fM->face_feature_idx[34][0] = 234;
			fM->face_feature_idx[35][0] = 7505;
			fM->face_feature_idx[36][0] = 5;
			fM->face_feature_idx[37][0] = 6844;
			fM->face_feature_idx[38][0] = 7136;

			//第一层--脸部
			myH->p_DeformHandler->DeformInit(v_id, V, V);
			myH->p_DeformHandler->setconset(facefixidx, facefixweights, facefixpos);
			myH->p_DeformHandler->addconset(fM->face_feature_idx, fM->face_feature_weights, fM->face_feature_pos);
			//myH->p_DeformHandler->setconset(layer1_idx, layer1_weights, layer1_pos);
			//myH->p_DeformHandler->addconset(layer2_idx, layer2_weights, layer2_pos);
			//myH->p_DeformHandler->addconset(layer3_idx, layer3_weights, layer3_pos);
			//myH->p_DeformHandler->DeformInit1(v_id, V, V);
			myH->p_DeformHandler->econ_threshold = 1.0;
			start = clock();
			myH->p_DeformHandler->Deform(v_id, V, V);
			finish = clock();
			totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
			timeRecord << totaltime << " ";
			myH->p_DeformHandler->ModVMat(V);
			
			sRV1 = V;
			sRV2 = s * (fromsmpltoPc * sRV1.transpose());
			sRV2.colwise() += trans;
			sRV3 = sRV2.transpose();
			sprintf_s(meshpath, 100, "data\\%d\\%d\\layer8.off", frame, i_th);
			fM->spCM->outMesh1(sRV3, meshpath);
			fM->spCM->outMesh1(sRV1, "refMesh2.off");

			
//			std::vector<std::vector<int>> faceidx1;
//			std::vector<std::vector<double>> faceweights1;
//			std::vector<Eigen::Vector3d> facepos1;
//			for (int i = 0; i < 20; i++)
//			{
//				faceidx1.push_back(fM->face_feature_idx[i]);
//				faceweights1.push_back(fM->face_feature_weights[i]);
//				facepos1.push_back(fM->face_feature_pos[i]);
//			}
//			std::vector<std::vector<int>> faceidx2;
//			std::vector<std::vector<double>> faceweights2;
//			std::vector<Eigen::Vector3d> facepos2;
//			for (int i = 20; i < 54; i++)
//			{
//				faceidx2.push_back(fM->face_feature_idx[i]);
//				faceweights2.push_back(fM->face_feature_weights[i]);
//				facepos2.push_back(fM->face_feature_pos[i]);
//			}
//			myH = new mns::MyMeshHandler("refMesh2.off", "parameters2.txt", argc, argv);
//			myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\male_expressions\\DE\\C.txt");
//			myH->p_DeformHandler->DeformInit(v_id, V, V);
//			myH->p_DeformHandler->setconset(facefixidx, facefixweights, facefixpos);
//			//myH->p_DeformHandler->addconset(faceidx2, faceweights2, facepos2);
//			myH->p_DeformHandler->addconset(faceidx1, faceweights1, facepos1);
//			myH->p_DeformHandler->econ_threshold = 1.0;
//			//myH->p_DeformHandler->Deform(v_id, V, V);
//			myH->p_DeformHandler->ModVMat(V);
//			sRV1 = V;
//			sRV2 = s * (fromsmpltoPc * sRV1.transpose());
//			sRV2.colwise() += trans;
//			sRV3 = sRV2.transpose();
//			/*sprintf_s(meshpath, 100, "data\\%d\\%d\\layer2.off", frame, i_th);
//			fM->spCM->outMesh1(sRV3, meshpath);*/
//
//			//continue;
//
//			//std::vector<int> layer5 = { 10387, 8808, 1239, 10822, 5136, 1640, 1325, 2998, 1780, 1642, 5081, 5397, 5606, 5661, 4992, 5391, 5037, 1586, 1615, 1549, 1969, 2221, 2287, 2051, 779, 947, 1102, 1075, 3294, 4814, 4818, 6665, 6681, 4492, 4559, 4549, 4588, 6526, 4904, 1457, 1049, 1464, 3154, 3304, 1004, 1094, 959, 1767};
//			//std::vector<std::vector<int>> layer5_idx(layer5.size());// +facecontour_pos.size());
//			//std::vector<std::vector<double>> layer5_weights(layer5.size());// +facecontour_pos.size());
//			//std::vector<Eigen::Vector3d> layer5_pos(layer5.size());// +facecontour_pos.size());
//
//			//for (int i = 0; i < layer5.size(); i++)
//			//{
//			//	vector<int> idx = { layer5[i] };
//			//	vector<double> weight = { 1.0 };
//			//	Eigen::Vector3d content = sRV1.row(layer5[i]).transpose();
//			//	layer5_idx[i] = idx;
//			//	layer5_weights[i] = weight;
//			//	layer5_pos[i] = content;
//			//}
//
//			////第二层
//			//myH->p_DeformHandler->setconset(layer5_idx, layer5_weights, layer5_pos);
//			//myH->p_DeformHandler->addconset(fM->face_feature_idx, fM->face_feature_weights, fM->face_feature_pos);
//			////myH->p_DeformHandler->addconset(layer4_idx, layer4_weights, layer4_pos);
//			//myH->p_DeformHandler->alpha_c = 800000;
//			//myH->p_DeformHandler->DeformInit1(v_id, V, V);
//			//myH->p_DeformHandler->econ_threshold = 0.2;
//			//myH->p_DeformHandler->Deform(v_id, V, V);
//			//myH->p_DeformHandler->ModVMat(V);
//			//sRV1 = V;
//			//sRV2 = s * (fromsmpltoPc * sRV1.transpose());
//			//sRV2.colwise() += trans;
//			//sRV3 = sRV2.transpose();
//			//sprintf_s(meshpath, 100, "data\\%d\\%d\\layer2.off", frame, i_th);
//			//fM->spCM->outMesh1(sRV3, meshpath);
//			//fM->spCM->setMesh(sRV1);
//			//fM->spCM->calculateFacesNormals();
//			//fM->spCM->calculateFacesFrame();
//			//fM->spCM->computeDiEdge();
//			//Eigen::VectorXd cm1 = fM->spCM->DiEdgeDataMatrix.row(0);
//			//Eigen::VectorXd cm2(cm1.size());
//			//for (int i = 0; i < fM->nE; i++)
//			//{
//			//	cm2(2 * i + 0) = cm1(2 * i + 0) - cm0(2 * i + 0);
//			//	cm2(2 * i + 1) = (cm1(2 * i + 1) / cm0(2 * i + 1)) - 1;
//			//}
//			//sprintf_s(meshpath, 100, "data\\%d\\%d\\cm.txt", frame, i_th);
//			//ofstream file_cm(meshpath);
//			//for (int i = 0; i < cm2.size(); i++)
//			//{
//			//	file_cm << cm2(i) << endl;
//			//}
//			//file_cm.close();
//			/*myH->p_DeformHandler->DeformInit1(v_id, V, V);
//			myH->p_DeformHandler->addconset(layer2_idx, layer2_weights, layer2_pos);
//			myH->p_DeformHandler->Deform(v_id, V, V);
//			myH->p_DeformHandler->ModVMat(V);
//			sRV = V.transpose();
//			sRV1 = sRV.transpose();
//			sRV2 = s * (fromsmpltoPc * sRV);
//			sRV2.colwise() += trans;
//			sRV3 = sRV2.transpose();
//			sprintf_s(meshpath, 100, "data\\%d\\%d\\layer2.off", frame, i_th);
//			fM->spCM->outMesh1(sRV3, meshpath);
//
//			myH->p_DeformHandler->DeformInit1(v_id, V, V);
//			myH->p_DeformHandler->addconset(layer3_idx, layer3_weights, layer3_pos);
//			myH->p_DeformHandler->Deform(v_id, V, V);
//			myH->p_DeformHandler->ModVMat(V);
//			sRV = V;
//			sRV1 = sRV.transpose();
//			sRV2 = s * (fromsmpltoPc * sRV1);
//			sRV2.colwise() += trans;
//			sRV3 = sRV2.transpose();
//			sprintf_s(meshpath, 100, "data\\%d\\%d\\layer3.off", frame, i_th);
//			fM->spCM->outMesh1(sRV3, meshpath);
//*/
//			/* 第二层 身体*/
//			PGMesh * pg0 = new PGMesh();
//			outfile1 = ofstream("initialize.txt");
//			for (int i = 0; i < 400; i++)
//			{
//				outfile1 << 0 << endl;
//			}
//			outfile1.close();
//			sprintf_s(meshpath, 100, "data\\%d\\%d\\layer2.off", frame, i_th);
//			OpenMesh::IO::read_mesh(*pg0, meshpath);
//			Eigen::VectorXd tempV = fM->getVertexMatrix(pg0);
//			delete pg0;
//			Eigen::MatrixXd V5(3, fM->nV);
//			for (int i = 0; i < fM->nV; i++)
//			{
//				V5.col(i) = tempV.segment(3 * i, 3);
//			}
//			V5.colwise() -= trans;
//			V5 = (1.0 / s) * V5;
//
//			/*sRV2 = fromPctosmpl * sRV2;
//			sRV = sRV2.transpose();*/
//			//after
//			Eigen::MatrixXd V6 = (fromPctosmpl * V5).transpose();
//			//sRV = sRV2.transpose();
//
//			fM->spCM->outMesh1(V6, "refMesh2.off");
//			Eigen::MatrixXd V7(fM->nV, 3);
//			myH = new mns::MyMeshHandler("refMesh2.off", "parameters1_DFAUST.txt", argc, argv);
//			//myH->p_DeformHandler->ReceiveBasis2("C.txt");
//			//myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\female_without_mouth\\DE\\C.txt");
//			myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\MFaust\\scripts\\male\\mixedDE\\C.txt");
//			myH->p_DeformHandler->DeformInit(v_id, V7, V7);
//			std::vector<int> bodyfix = {3726, 9598, 3591, 9090, 3592, 9088, 9865, 9866, 9110, 9108, 9312, 10706, 9309, 9859, 10048, 10044, 9953, 3609, 9110, 3593, 9090, 9092, 7032, 7882, 7974, 7029, 8601, 7025, 351, 7866, 8624, 8610, 7616, 2786, 8622, 8648, 8633,7207, 25, 9768, 9047, 10614, 8755, 6973, 86, 7242, 6911, 8788, 9438, 8990, 9041,9272, 9714, 9883, 10602, 8766, 7059, 7294, 139, 268, 7537, 7834, 7838, 9959, 9681, 3869, 9334, 9665, 9818, 9954, 9546, 9531, 9775, 7646, 7653, 245, 159, 7460, 7012, 129, 375, 163, };
//			std::vector<std::vector<int>> bodyfixidx;
//			std::vector<std::vector<double>> bodyfixweights;
//			std::vector<Eigen::Vector3d> bodyfixpos;
//			for (int i = 0; i < bodyfix.size(); i++)
//			{
//				vector<int> idx = { bodyfix[i] };
//				vector<double> weight = { 1.0 };
//				Eigen::Vector3d coord;
//				coord = V6.row(bodyfix[i]).transpose();
//				bodyfixidx.push_back(idx);
//				bodyfixpos.push_back(coord);
//				bodyfixweights.push_back(weight);
//			}
//			myH->p_DeformHandler->setconset(bodyfixidx, bodyfixweights, bodyfixpos);
//			myH->p_DeformHandler->addconset(layer1_idx, layer1_weights, layer1_pos);
//			//myH->p_DeformHandler->addconset(layer2_idx, layer2_weights, layer2_pos);
//			//myH->p_DeformHandler->DeformInit1(v_id, V, V);
//			myH->p_DeformHandler->econ_threshold = 1.5;
//			start = clock();
//			//myH->p_DeformHandler->Deform(v_id, V7, V7);
//			finish = clock();
//			totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
//			timeRecord << totaltime << " ";
//			myH->p_DeformHandler->ModVMat(V7);
//			
//			Eigen::MatrixXd V8 = s * (fromsmpltoPc * V7.transpose());
//			V8.colwise() += trans;
//			Eigen::MatrixXd V9 = V8.transpose();
//			/*sprintf_s(meshpath, 100, "data\\%d\\%d\\layer3.off", frame, i_th);
//			fM->spCM->outMesh1(V9, meshpath);*/
//			
//			std::vector<int> bodyfix4 = { 3049, 4212, 1460, 4904, 4365, 877, 628, 1394, 5151, 4089, 7591, 1463, 1168, 4496, 4481, 1675, 1314, 4768, 5117 };
//			std::vector<std::vector<int>> bodyfixidx4;
//			std::vector<std::vector<double>> bodyfixweights4;
//			std::vector<Eigen::Vector3d> bodyfixpos4;
//			for (int i = 0; i < bodyfix4.size(); i++)
//			{
//				vector<int> idx = { bodyfix4[i] };
//				vector<double> weight = { 1.0 };
//				Eigen::Vector3d coord;
//				coord = V7.row(bodyfix4[i]).transpose();
//				bodyfixidx4.push_back(idx);
//				bodyfixpos4.push_back(coord);
//				bodyfixweights4.push_back(weight);
//			}
//			myH->p_DeformHandler->setconset(bodyfixidx, bodyfixweights, bodyfixpos);
//			myH->p_DeformHandler->addconset(layer1_idx, layer1_weights, layer1_pos);
//			myH->p_DeformHandler->addconset(layer2_idx, layer2_weights, layer2_pos);
//			myH->p_DeformHandler->DeformInit1(v_id, V7, V7);
//			//myH->p_DeformHandler->Deform(v_id, V7, V7);
//			myH->p_DeformHandler->ModVMat(V7);
//
//			V8 = s * (fromsmpltoPc * V7.transpose());
//			V8.colwise() += trans;
//			V9 = V8.transpose();
//			/*sprintf_s(meshpath, 100, "data\\%d\\%d\\layer4.off", frame, i_th);
//			fM->spCM->outMesh1(V9, meshpath);*/


			//手部姿态
			PGMesh * pg1 = new PGMesh();
			outfile1 = ofstream("initialize.txt");
			for (int i = 0; i < 100; i++)
			{
				outfile1 << 0 << endl;
			}
			outfile1.close();
			sprintf_s(meshpath, 100, "data\\%d\\%d\\layer4.off", frame, i_th);
			OpenMesh::IO::read_mesh(*pg1, meshpath);
			Eigen::VectorXd tempV1 = fM->getVertexMatrix(pg1);
			delete pg1;
			Eigen::MatrixXd V10(3, fM->nV);
			for (int i = 0; i < fM->nV; i++)
			{
				V10.col(i) = tempV1.segment(3 * i, 3);
			}
			V10.colwise() -= trans;
			V10 = (1.0 / s) * V10;
			/*sRV2 = fromPctosmpl * sRV2;
			sRV = sRV2.transpose();*/
			//after
			Eigen::MatrixXd V11 = (fromPctosmpl * V10).transpose();
			//Eigen::MatrixXd V11 = V10.transpose();
			//sRV = sRV2.transpose();
			fM->spCM->outMesh1(V11, "refmesh.off");
			fM->spCM->outMesh1(V11, "refMesh2.off");
			Eigen::MatrixXd V12(fM->nV, 3);
			myH = new mns::MyMeshHandler("refMesh2.off", "parameters3.txt", argc, argv);
			myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\male_hand\\DE\\newC.txt");
			myH->p_DeformHandler->DeformInit(v_id, V12, V12);
			std::vector<int> handfix = { 4087, 5177, 3038, 2896, 3523, 1678, 1668, 10450, 5221, 4840, 1709 };
			std::vector<std::vector<int>> handfixidx;
			std::vector<std::vector<double>> handfixweights;
			std::vector<Eigen::Vector3d> handfixpos;
			for (int i = 0; i < handfix.size(); i++)
			{
				vector<int> idx = { handfix[i] };
				vector<double> weight = { 1.0 };
				Eigen::Vector3d coord;
				coord = V11.row(handfix[i]).transpose();
				handfixidx.push_back(idx);
				handfixpos.push_back(coord);
				handfixweights.push_back(weight);
			}
			myH->p_DeformHandler->setconset(handfixidx, handfixweights, handfixpos);
			/*myH->p_DeformHandler->setconset(layer1_idx, layer1_weights, layer1_pos);
			myH->p_DeformHandler->addconset(layer2_idx, layer2_weights, layer2_pos);
			myH->p_DeformHandler->addconset(layer3_idx, layer3_weights, layer3_pos);*/
			myH->p_DeformHandler->addconset(fM->hand_feature_idx, fM->hand_feature_weights, fM->hand_feature_pos);
			myH->p_DeformHandler->econ_threshold = 3;
			//myH->p_DeformHandler->DeformInit1(v_id, V, V);
			start = clock();
			myH->p_DeformHandler->Deform(v_id, V12, V12);
			//myH->p_DeformHandler->Deform(v_id, V12, V12);
			finish = clock();
			totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
			timeRecord << totaltime << endl;
			myH->p_DeformHandler->ModVMat(V12);

			Eigen::MatrixXd V13 = s * (fromsmpltoPc * V12.transpose());
			V13.colwise() += trans;
			Eigen::MatrixXd V14 = V13.transpose();
			sprintf_s(meshpath, 100, "data\\%d\\%d\\layer5.off", frame, i_th);
			fM->spCM->outMesh1(V14, meshpath);
			sprintf_s(meshpath, 100, "data\\%d.off", frame);
			fM->spCM->outMesh1(V14, meshpath);


			//手部姿态2
			PGMesh * pg2 = new PGMesh();
			outfile1 = ofstream("initialize.txt");
			for (int i = 0; i < 100; i++)
			{
				outfile1 << 0 << endl;
			}
			outfile1.close();
			sprintf_s(meshpath, 100, "data\\%d\\%d\\layer5.off", frame, i_th);
			OpenMesh::IO::read_mesh(*pg2, meshpath);
			Eigen::VectorXd tempV2 = fM->getVertexMatrix(pg2);
			delete pg2;
			Eigen::MatrixXd V15(3, fM->nV);
			for (int i = 0; i < fM->nV; i++)
			{
				V15.col(i) = tempV2.segment(3 * i, 3);
			}
			V15.colwise() -= trans;
			V15 = (1.0 / s) * V15;
			/*sRV2 = fromPctosmpl * sRV2;
			sRV = sRV2.transpose();*/
			//after
			Eigen::MatrixXd V16 = (fromPctosmpl * V15).transpose();
			//Eigen::MatrixXd V11 = V10.transpose();
			//sRV = sRV2.transpose();
			fM->spCM->outMesh1(V16, "refmesh.off");
			fM->spCM->outMesh1(V16, "refMesh2.off");
			Eigen::MatrixXd V17(fM->nV, 3);
			myH = new mns::MyMeshHandler("refMesh2.off", "parameters3.txt", argc, argv);
			myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\male_hand\\DE\\newC.txt");
			myH->p_DeformHandler->DeformInit(v_id, V17, V17);
			
			
			myH->p_DeformHandler->setconset(handfixidx, handfixweights, handfixpos);
			/*myH->p_DeformHandler->setconset(layer1_idx, layer1_weights, layer1_pos);
			myH->p_DeformHandler->addconset(layer2_idx, layer2_weights, layer2_pos);
			myH->p_DeformHandler->addconset(layer3_idx, layer3_weights, layer3_pos);*/
			myH->p_DeformHandler->addconset(fM->hand_feature_idx, fM->hand_feature_weights, fM->hand_feature_pos);
			myH->p_DeformHandler->econ_threshold = 3;
			//myH->p_DeformHandler->DeformInit1(v_id, V, V);
			start = clock();
			myH->p_DeformHandler->Deform(v_id, V17, V17);
			//myH->p_DeformHandler->Deform(v_id, V12, V12);
			finish = clock();
			totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
			timeRecord << totaltime << endl;
			myH->p_DeformHandler->ModVMat(V17);

			Eigen::MatrixXd V18 = s * (fromsmpltoPc * V17.transpose());
			V18.colwise() += trans;
			Eigen::MatrixXd V19 = V18.transpose();
			sprintf_s(meshpath, 100, "data\\%d\\%d\\layer6.off", frame, i_th);
			fM->spCM->outMesh1(V19, meshpath);
			sprintf_s(meshpath, 100, "data\\%d.off", frame);
			fM->spCM->outMesh1(V19, meshpath);

		}
	}

	//cout << y;
}

void testReconstSeq(int argc, char** argv)
{
	clock_t start, finish;
	char *meshpath = new char[100];
	sprintf_s(meshpath, 100, "F:\\scan\\male_without_mouth\\0000.obj");
	PGMesh *pg = new PGMesh();
	OpenMesh::IO::read_mesh(*pg, meshpath);
	Eigen::VectorXd y(69);
	Eigen::VectorXd x(10);
	Eigen::VectorXd facecoef(200);
	Eigen::VectorXd gesturecoef(100);
	Eigen::VectorXd poseparm(500);
	Eigen::Matrix3d globalrotate, fromPctosmpl, fromsmpltoPc;
	double s;
	std::vector<int> layer = { 0, 1, 2, 3, 4, 5,6,7,8,9,10,11,12,13,14,15,16,17,18,19, 20, 21, 22, 23, 24, 25, 26, 27, 28 };
	std::vector<int> transvec = { 0, 6, 7, 8, 12, 13, 14 };
	std::vector<int> transvec1 = { 5, 8, 11, 14 };
	//std::vector<int> layer1 = { 3, 6, 9, 12 , 19, 20, 21, 22, 23, 24, 25, 26 };
	std::vector<int> layer1 = { 3, 6, 9, 12  };
	std::vector<int> layer2 = { 4, 7, 10, 13, 5, 8, 11, 14 };
	std::vector<int> layer3 = { 5, 8, 11, 14 };
	Eigen::Vector3d  trans, global_trans;
	FaceBodyModel * fM = new FaceBodyModel(10, 250, pg);
	Eigen::VectorXd out(fM->nV * 3);
	sprintf_s(meshpath, 100, "F:\\scan\\male_expressions\\DE\\newC.txt");
	fM->loadExpressionBases(meshpath, 200);
	sprintf_s(meshpath, 100, "F:\\project\\splocs-master\\newC.txt");
	//sprintf_s(meshpath, 100, "F:\\scan\\female_without_mouth\\DE\\C.txt");
	//sprintf_s(meshpath, 100, "E:\\reconstruction\\250\\work\\C.txt");
	fM->loadPoseBases(meshpath, 400);
	sprintf_s(meshpath, 100, "F:\\scan\\male_hand\\DE\\newC.txt");
	fM->loadgestureBases(meshpath, 100);
	sprintf_s(meshpath, 100, "F:\\scan\\smpl\\male_shape\\mean.txt");
	fM->loadMeanShape(meshpath);
	sprintf_s(meshpath, 100, "F:\\scan\\smpl\\male_shape\\shapebasis.txt");
	fM->loadShapeBases(meshpath);

	//output shape basis
	//for (int i = 0; i < 10; i++)
	//{
	//	Eigen::VectorXd s_coef(10);
	//	s_coef.setZero();
	//	s_coef(i) = 3.0;
	//	Eigen::VectorXd y_out = fM->generateShapeWithoutPose(s_coef);
	//	std::vector<double> error(fM->nV);
	//	for (int j = 0; j < fM->nV; j++)
	//	{
	//		error[j] = (fM->C_S.block<1, 3>(i, 3 * j + 0)) * (fM->C_S.block<1, 3>(i, 3 * j + 0).transpose());
	//		error[j] = sqrt(error[j]);
	//	}
	//	sprintf_s(meshpath, 100, "shape_component\\mesh_%04d.obj", i);
	//	fM->spCM->colorize(error, y_out, meshpath);
	//}

	Eigen::MatrixXd deformedV(fM->nV, 3);
	fM->readBodyFeatureidx("body_feature_idx.txt");
	int a[19] = { 12, 15, 0, 16, 18, 20, 1, 4, 7, 17, 19, 21, 2, 5, 8, 24, 25, 26, 27 };
	fM->selectbodyfeature.resize(19);
	copy(a, a + 19, fM->selectbodyfeature.begin());
	fM->body_feature_idx1.resize(19);
	fM->body_feature_weights1.resize(19);
	/*fM->body_feature_idx1.resize(29);
	fM->body_feature_weights1.resize(29);*/
	for (int i = 0; i < 19; i++)
	{
		fM->body_feature_idx1[i] = fM->body_feature_idx[a[i]];
		fM->body_feature_weights1[i] = fM->body_feature_weights[a[i]];
	}
	/*std::vector<int> bodyadd_conset_idx = { 570, 8359, 4152, 5207, 751, 2994, 4093, 677, 4097, 1861 };
	for (int i = 19; i < 29; i++)
	{
		std::vector<int> idx = { bodyadd_conset_idx[i - 19] };
		std::vector<double> weight = { 1.0 };
	}*/
	std::auto_ptr<mns::MyMeshHandler> p_meshHandler(new mns::MyMeshHandler(meshpath, fM->body_feature_idx1, fM->body_feature_weights1));
	Eigen::VectorXd tmpV(fM->nV * 3);
	std::vector<int> smplbody0 = { 3500, 4985, 912, 5333, 1294, 3022, 2967, 6427, 3107, 6530,1619, 1665, 5212, 5165, 4515, 1029, 4486, 1001 };
	std::vector<int> fbmmbody0 = { 4785, 4956, 912, 5306, 1294, 2995, 2940, 6373, 3080, 6476 ,1619, 1665, 5185, 5138, 4488, 1029, 4590, 1001 };
	std::vector<int> smplbody1 = { 3500, 4985, 912, 5333, 1294, 3022, 2967, 6427, 3107, 6530 };
	std::vector<int> fbmmbody1 = { 4785, 4956, 912, 5306, 1294, 2995, 2940, 6373, 3080, 6476 };
	std::vector<int> smplbody2 = { 1619, 1665, 5212, 5165, 4515, 1029, 4486, 1001 };
	std::vector<int> fbmmbody2 = { 1619, 1665, 5185, 5138, 4488, 1029, 4590, 1001 };
	std::vector<int> smplbody3 = { 1970, 1930, 5669, 5560, 6723, 6581, 3323, 3203 };
	std::vector<int> fbmmbody3 = { 1970, 1930, 5642, 5533, 6669, 6527, 3296, 3176 };
	PGMesh * tmpmesh = new PGMesh();
	/*vector<int> frameset(8499-8399+1);
	int start_frame = 8399;
	for (int i = 0; i < 8499 - 8399 + 1; i++)
	{
		frameset[i] = 8399 + i;
	}*/
	vector<int> frameset = { 8372};
	
	ofstream timeRecord("time_statistic.txt");
	for (int ss = 0; ss < frameset.size(); ss++)
	{
		int frame = frameset[ss];
		vector<int> ith_set = { 2 };

		//for (int sss = 0; sss < ith_set.size(); sss++)
		//{
		//	int i_th = ith_set[sss];
		//	sprintf_s(meshpath, 100, "data\\%d\\%d\\layer5.off", frame, i_th);
		//	PGMesh *pg000 = new PGMesh();
		//	OpenMesh::IO::read_mesh(*pg000, meshpath);
		//	sprintf_s(meshpath, 100, "data\\mesh%d_%d.obj", frame, i_th);
		//	OpenMesh::IO::write_mesh(*pg000, meshpath);
		//}
		
		for (int sss = 0; sss < ith_set.size(); sss++)
		{
			for (int newstate = 0; newstate < 1; newstate++)
			{
				int i_th = ith_set[sss];
				/*sprintf_s(meshpath, 100, "data\\%d\\%d\\layer5.off", frame, i_th);
				PGMesh *pg000 = new PGMesh();
				OpenMesh::IO::read_mesh(*pg000, meshpath);*/
				/*sprintf_s(meshpath, 100, "data\\mesh%d_%d.obj", frame, i_th);
				OpenMesh::IO::write_mesh(*pg000, meshpath);
				continue;*/
				sprintf_s(meshpath, 100, "data\\%d\\%d\\coef.txt", frame, i_th);
				readcoefficient(meshpath, x, y);
				sprintf_s(meshpath, 100, "data\\%d\\%d\\global.txt", frame, i_th);
				read_global(meshpath, globalrotate, fromPctosmpl, fromsmpltoPc, s, trans, global_trans);
				sprintf_s(meshpath, 100, "data\\%d\\%d\\anchorpoint.txt", frame, i_th);
				fM->spCM->loadAnchorData(meshpath);
				fM->spCM->presolve();
				sprintf_s(meshpath, 100, "data\\%d\\%d\\bodyconset.txt", frame, i_th);
				//fM->body_feature_pos1.resize(29);
				fM->body_feature_pos1.resize(19);
				readBodyFeaturePos(meshpath, fM->body_feature_pos1);
				fM->body_feature_pos.resize(19);
				//fM->body_feature_pos.resize(29);
				for (int i = 0; i < 19; i++)
				{
					fM->body_feature_pos[i] = fM->body_feature_pos1[i];
				}
				sprintf_s(meshpath, 100, "data\\%d\\%d\\handconset.txt", frame, i_th);
				fM->hand_feature_idx.resize(42);
				fM->hand_feature_weights.resize(42);
				fM->hand_feature_pos.resize(42);
				readHandFeaturePos(meshpath, fM->hand_feature_pos, fM->hand_feature_idx, fM->hand_feature_weights);
				sprintf_s(meshpath, 100, "data\\%d\\%d\\faceconset.txt", frame, i_th);
				fM->face_feature_idx.resize(54);
				fM->face_feature_pos.resize(54);
				fM->face_feature_weights.resize(54);
				readFaceFeaturePos(meshpath, fM->face_feature_pos, fM->face_feature_idx, fM->face_feature_weights);
				int b[54] = { 9160, 3030, 152, 7173, 9811, 9249, 7258, 8782, 444, 4250, 5242, 1871, 1574, 6445, 4778, 5237, 932, 4336, 4407, 1160, 3098, 1099, 4613, 4578, 3294, 1127, 1050, 870, 4313, 4510, 6551, 3131, 3146, 1216, 2952, 1909, 1976, 6303, 5141, 5401, 5499, 5619, 2794, 1569, 7629, 7690, 4759, 1244, 4167, 6256, 1435, 890, 3063, 4303 };
				/*for (int i = 0; i < 54; i++)
				{
				fM->face_feature_idx[54 + i].resize(1);
				fM->face_feature_idx[54 + i][0] = b[i];
				fM->face_feature_idx[54 + i][0] = b[i];
				fM->face_feature_weights[54 + i].resize(1);
				fM->face_feature_weights[54 + i][0] = 1.0;
				}*/

				/*std::vector<std::vector<int>> layer1_idx(fbmmbody1.size());
				std::vector<std::vector<double>> layer1_weights(fbmmbody1.size());
				std::vector<Eigen::Vector3d> layer1_pos(fbmmbody1.size());

				std::vector<std::vector<int>> layer2_idx(fbmmbody2.size() + 4);
				std::vector<std::vector<double>> layer2_weights(fbmmbody2.size()+ 4);
				std::vector<Eigen::Vector3d> layer2_pos(fbmmbody2.size()+ 4);

				std::vector<std::vector<int>> layer3_idx(fbmmbody3.size());
				std::vector<std::vector<double>> layer3_weights(fbmmbody3.size());
				std::vector<Eigen::Vector3d> layer3_pos(fbmmbody3.size());*/

				std::vector<std::vector<int>> layer1_idx(layer1.size());
				std::vector<std::vector<double>> layer1_weights(layer1.size());
				std::vector<Eigen::Vector3d> layer1_pos(layer1.size());

				std::vector<std::vector<int>> layer2_idx(layer2.size());
				std::vector<std::vector<double>> layer2_weights(layer2.size());
				std::vector<Eigen::Vector3d> layer2_pos(layer2.size());

				std::vector<std::vector<int>> layer3_idx(layer3.size());
				std::vector<std::vector<double>> layer3_weights(layer3.size());
				std::vector<Eigen::Vector3d> layer3_pos(layer3.size());

				std::vector<std::vector<int>> facecontour_idx;
				std::vector<std::vector<double>> facecontour_weights;
				std::vector<Eigen::Vector3d> facecontour_pos;
				sprintf_s(meshpath, 100, "data\\%d\\%d\\facecontourconset.txt", frame, i_th);
				readFacecontourFeaturePos(meshpath, facecontour_pos, facecontour_idx, facecontour_weights);

				std::vector<int> layer4 = { 15, 16, 17, 18 };
				std::vector<std::vector<int>> layer4_idx(layer4.size());// +facecontour_pos.size());
				std::vector<std::vector<double>> layer4_weights(layer4.size());// +facecontour_pos.size());
				std::vector<Eigen::Vector3d> layer4_pos(layer4.size());// +facecontour_pos.size());

																	   /*for (int i = 0; i < facecontour_idx.size(); i++)
																	   {
																	   fM->face_feature_idx.push_back(facecontour_idx[i]);
																	   fM->face_feature_weights.push_back(facecontour_weights[i]);
																	   fM->face_feature_pos.push_back(facecontour_pos[i]);
																	   }*/
				start = clock();
				poseparm = fM->getPoseParm(x, y);
				finish = clock();
				double totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
				timeRecord << totaltime << " ";
				ofstream outfile1("initialize.txt");
				for (int i = 0; i < 200; i++)
				{
					outfile1 << 0 << endl;
				}
				for (int i = 0; i < 200; i++)
				{
					outfile1 << 0 << endl;
				}
				outfile1.close();
				//cout << poseparm << endl;
				start = clock();
				Eigen::MatrixXd V = fM->f_x(x, poseparm);
				finish = clock();
				totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
				timeRecord << totaltime << " ";
				Eigen::VectorXd ref = fM->generateShapeWithoutPose(x);
				Eigen::MatrixXd refMatrix(fM->nV, 3);
				for (int i = 0; i < fM->nV; i++)
				{
					refMatrix.row(i) = ref.segment(3 * i, 3);
				}
				//fM->spCM->outMesh1(ref, "refmesh.off");
				sprintf_s(meshpath, 100, "data\\%d\\%d\\deformed.obj", frame, i_th);
				OpenMesh::IO::read_mesh(*tmpmesh, meshpath);
				Eigen::MatrixXd smplposMatrix = fM->getVertexMatrix1(tmpmesh);
				smplposMatrix.rowwise() -= trans.transpose();
				Eigen::MatrixXd smplposMatrix_T = (1.0 / s) * (smplposMatrix.transpose());
				smplposMatrix_T = fromPctosmpl * (smplposMatrix_T);
				smplposMatrix = (smplposMatrix_T.transpose());

				Eigen::MatrixXd smplFeaturepos(smplbody0.size(), 3);
				for (int i = 0; i < smplbody0.size(); i++)
				{
					smplFeaturepos.row(i) = smplposMatrix.row(smplbody0[i]);
				}
				Eigen::MatrixXd smplFeaturepos1(smplbody1.size(), 3);
				Eigen::MatrixXd smplFeaturepos2(smplbody2.size(), 3);
				Eigen::MatrixXd smplFeaturepos3(smplbody3.size(), 3);
				std::vector<std::vector<int>> fbmmbody1_idx(smplbody1.size());
				std::vector<std::vector<double>> fbmmbody1_weight(smplbody1.size());
				std::vector<std::vector<int>> fbmmbody2_idx(smplbody2.size());
				std::vector<std::vector<double>> fbmmbody2_weight(smplbody2.size());
				std::vector<std::vector<int>> fbmmbody3_idx(smplbody3.size());
				std::vector<std::vector<double>> fbmmbody3_weight(smplbody3.size());
				for (int i = 0; i < smplbody1.size(); i++)
				{
					smplFeaturepos1.row(i) = smplposMatrix.row(smplbody1[i]);
					std::vector<int> idx = { smplbody1[i] };
					std::vector<double> weight = { 1.0 };
					fbmmbody1_idx[i] = idx;
					fbmmbody1_weight[i] = weight;
				}
				for (int i = 0; i < smplbody2.size(); i++)
				{
					smplFeaturepos2.row(i) = smplposMatrix.row(smplbody2[i]);
					std::vector<int> idx = { smplbody2[i] };
					std::vector<double> weight = { 1.0 };
					fbmmbody2_idx[i] = idx;
					fbmmbody2_weight[i] = weight;
				}
				for (int i = 0; i < smplbody3.size(); i++)
				{
					smplFeaturepos3.row(i) = smplposMatrix.row(smplbody3[i]);
					std::vector<int> idx = { smplbody1[3] };
					std::vector<double> weight = { 1.0 };
					fbmmbody3_idx[i] = idx;
					fbmmbody3_weight[i] = weight;
				}
				//Eigen::VectorXd currentFeaturepos = fM->getBodyPosfromCurrentV1(V, fbmmbody0);
				Eigen::VectorXd currentFeaturepos = fM->getFacePosfromCurrentV(V);
				//优化全局旋转参数R,t
				//fM->findRigidTransform(currentFeaturepos, globalrotate, global_trans, smplFeaturepos);
				fM->findRigidTransform(currentFeaturepos, globalrotate, global_trans, fM->face_feature_pos);
				//Eigen::VectorXd currentFeaturepos = fM->getBodyPosfromCurrentV(V, layer);
				////			//优化全局旋转参数R,t
				//fM->findRigidTransform(currentFeaturepos, globalrotate, global_trans, fM->body_feature_pos);
				refMatrix = (globalrotate * refMatrix.transpose()).transpose();
				refMatrix.rowwise() += global_trans.transpose();
				fM->spCM->outMesh1(refMatrix, "refmesh.off");
				sprintf_s(meshpath, 100, "data\\%d\\%d\\ref.off", frame, i_th);
				fM->spCM->outMesh1(refMatrix, meshpath);
				fM->spCM->setMesh(refMatrix);
				fM->spCM->calculateFacesNormals();
				fM->spCM->calculateFacesFrame();
				fM->spCM->computeDiEdge();
				Eigen::VectorXd cm0 = fM->spCM->DiEdgeDataMatrix.row(0);

				for (int i = 0; i < layer1.size(); i++)
				{
					std::vector<int> idx = fM->body_feature_idx1[layer1[i]];
					std::vector<double> weight = fM->body_feature_weights1[layer1[i]];
					Eigen::Vector3d pos = fM->body_feature_pos[layer1[i]];
					layer1_idx[i] = idx;
					layer1_weights[i] = weight;
					layer1_pos[i] = pos;
				}
				for (int i = 0; i < layer2.size(); i++)
				{
					std::vector<int> idx = fM->body_feature_idx1[layer2[i]];
					std::vector<double> weight = fM->body_feature_weights1[layer2[i]];
					Eigen::Vector3d pos = fM->body_feature_pos[layer2[i]];
					layer2_idx[i] = idx;
					layer2_weights[i] = weight;
					layer2_pos[i] = pos;
				}
				for (int i = 0; i < layer3.size(); i++)
				{
					std::vector<int> idx = fM->body_feature_idx1[layer3[i]];
					std::vector<double> weight = fM->body_feature_weights1[layer3[i]];
					Eigen::Vector3d pos = fM->body_feature_pos[layer3[i]];
					layer3_idx[i] = idx;
					layer3_weights[i] = weight;
					layer3_pos[i] = pos;
				}

				for (int i = 0; i < layer4.size(); i++)
				{
					std::vector<int> idx = fM->body_feature_idx1[layer4[i]];
					std::vector<double> weight = fM->body_feature_weights1[layer4[i]];
					Eigen::Vector3d pos = fM->body_feature_pos[layer4[i]];
					layer4_idx[i] = idx;
					layer4_weights[i] = weight;
					layer4_pos[i] = pos;
				}


				Eigen::MatrixXd sRV = globalrotate * (V.transpose());
				sRV.colwise() += global_trans;
				Eigen::MatrixXd sRV1 = sRV.transpose();
				//currentFeaturepos = fM->getBodyPosfromCurrentV(sRV1, transvec);

				Eigen::MatrixXd sRV2 = s * (fromsmpltoPc * sRV1.transpose());
				sRV2.colwise() += trans;
				Eigen::MatrixXd sRV3 = sRV2.transpose();
				sprintf_s(meshpath, 100, "data\\%d\\%d\\initial_global.off", frame, i_th);
				fM->spCM->outMesh1(sRV3, meshpath);
				sprintf_s(meshpath, 100, "data\\%d\\%d\\initial_reconstruction.off", frame, i_th);
				fM->spCM->outMesh1(sRV1, meshpath);
				fM->spCM->outMesh1(sRV1, "refMesh2.off");
				p_meshHandler->fitmode = 0;

				mns::MyMeshHandler* myH = new mns::MyMeshHandler("refMesh2.off", "parameters2.txt", argc, argv);
				myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\male_expressions\\DE\\C.txt");
				//myH->p_DeformHandler->ReceiveBasis2("C.txt");
				std::vector<integer> v_id;

				std::vector<int> facefix = { 5298, 1864, 6448, 8387, 4782, 5105, 1651 };
				std::vector<std::vector<int>> facefixidx;
				std::vector<std::vector<double>> facefixweights;
				std::vector<Eigen::Vector3d> facefixpos;
				for (int i = 0; i < facefix.size(); i++)
				{
					vector<int> idx = { facefix[i] };
					vector<double> weight = { 1.0 };
					Eigen::Vector3d coord;
					coord = sRV1.row(facefix[i]).transpose();
					facefixidx.push_back(idx);
					facefixpos.push_back(coord);
					facefixweights.push_back(weight);
				}
				//第一层--脸部
				myH->p_DeformHandler->DeformInit(v_id, V, V);
				myH->p_DeformHandler->setconset(facefixidx, facefixweights, facefixpos);
				myH->p_DeformHandler->addconset(fM->face_feature_idx, fM->face_feature_weights, fM->face_feature_pos);
				//myH->p_DeformHandler->setconset(layer1_idx, layer1_weights, layer1_pos);
				//myH->p_DeformHandler->addconset(layer2_idx, layer2_weights, layer2_pos);
				//myH->p_DeformHandler->addconset(layer3_idx, layer3_weights, layer3_pos);
				//myH->p_DeformHandler->DeformInit1(v_id, V, V);
				myH->p_DeformHandler->econ_threshold = 1.0;
				start = clock();
				//if (newstate == 1)
				{
					myH->p_DeformHandler->Deform(v_id, V, V);
				}
				finish = clock();
				totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
				timeRecord << totaltime << " ";
				myH->p_DeformHandler->ModVMat(V);

				sRV1 = V;
				sRV2 = s * (fromsmpltoPc * sRV1.transpose());
				sRV2.colwise() += trans;
				sRV3 = sRV2.transpose();
				sprintf_s(meshpath, 100, "data\\%d\\%d\\layer1.off", frame, i_th);
				fM->spCM->outMesh1(sRV3, meshpath);
				fM->spCM->outMesh1(sRV1, "refMesh2.off");


				std::vector<std::vector<int>> faceidx1;
				std::vector<std::vector<double>> faceweights1;
				std::vector<Eigen::Vector3d> facepos1;
				for (int i = 0; i < 20; i++)
				{
					faceidx1.push_back(fM->face_feature_idx[i]);
					faceweights1.push_back(fM->face_feature_weights[i]);
					facepos1.push_back(fM->face_feature_pos[i]);
				}
				std::vector<std::vector<int>> faceidx2;
				std::vector<std::vector<double>> faceweights2;
				std::vector<Eigen::Vector3d> facepos2;
				for (int i = 20; i < 54; i++)
				{
					faceidx2.push_back(fM->face_feature_idx[i]);
					faceweights2.push_back(fM->face_feature_weights[i]);
					facepos2.push_back(fM->face_feature_pos[i]);
				}
				myH = new mns::MyMeshHandler("refMesh2.off", "parameters2.txt", argc, argv);
				myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\male_expressions\\DE\\C.txt");
				myH->p_DeformHandler->DeformInit(v_id, V, V);
				myH->p_DeformHandler->setconset(facefixidx, facefixweights, facefixpos);
				//myH->p_DeformHandler->addconset(faceidx2, faceweights2, facepos2);
				myH->p_DeformHandler->addconset(faceidx1, faceweights1, facepos1);
				myH->p_DeformHandler->econ_threshold = 1.0;
				//myH->p_DeformHandler->Deform(v_id, V, V);
				myH->p_DeformHandler->ModVMat(V);
				sRV1 = V;
				sRV2 = s * (fromsmpltoPc * sRV1.transpose());
				sRV2.colwise() += trans;
				sRV3 = sRV2.transpose();
				sprintf_s(meshpath, 100, "data\\%d\\%d\\layer2.off", frame, i_th);
				fM->spCM->outMesh1(sRV3, meshpath);


				/* 第二层 身体*/
				PGMesh * pg0 = new PGMesh();
				outfile1 = ofstream("initialize.txt");
				for (int i = 0; i < 500; i++)
				{
					outfile1 << 0 << endl;
				}
				outfile1.close();
				sprintf_s(meshpath, 100, "data\\%d\\%d\\layer2.off", frame, i_th);
				OpenMesh::IO::read_mesh(*pg0, meshpath);
				Eigen::VectorXd tempV = fM->getVertexMatrix(pg0);
				delete pg0;
				Eigen::MatrixXd V5(3, fM->nV);
				for (int i = 0; i < fM->nV; i++)
				{
					V5.col(i) = tempV.segment(3 * i, 3);
				}
				V5.colwise() -= trans;
				V5 = (1.0 / s) * V5;

				/*sRV2 = fromPctosmpl * sRV2;
				sRV = sRV2.transpose();*/
				//after
				Eigen::MatrixXd V6 = (fromPctosmpl * V5).transpose();
				//sRV = sRV2.transpose();

				fM->spCM->outMesh1(V6, "refMesh2.off");
				Eigen::MatrixXd V7(fM->nV, 3);
				myH = new mns::MyMeshHandler("refMesh2.off", "parameters1.txt", argc, argv);
				//myH->p_DeformHandler->ReceiveBasis2("C.txt");
				//myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\female_without_mouth\\DE\\C.txt");
				myH->p_DeformHandler->ReceiveBasis2("E:\\reconstruction\\250\\work\\newC.txt");
				myH->p_DeformHandler->DeformInit(v_id, V7, V7);
				std::vector<int> bodyfix = { 3726, 9598, 3591, 9090, 3592, 9088, 9865, 9866, 9110, 9108, 9312, 10706, 9309, 9859, 10048, 10044, 9953, 3609, 9110, 3593, 9090, 9092, 7032, 7882, 7974, 7029, 8601, 7025, 351, 7866, 8624, 8610, 7616, 2786, 8622, 8648, 8633,7207, 25, 9768, 9047, 10614, 8755, 6973, 86, 7242, 6911, 8788, 9438, 8990, 9041,9272, 9714, 9883, 10602, 8766, 7059, 7294, 139, 268, 7537, 7834, 7838, 9959, 9681, 3869, 9334, 9665, 9818, 9954, 9546, 9531, 9775, 7646, 7653, 245, 159, 7460, 7012, 129, 375, 163, };
				std::vector<std::vector<int>> bodyfixidx;
				std::vector<std::vector<double>> bodyfixweights;
				std::vector<Eigen::Vector3d> bodyfixpos;
				for (int i = 0; i < bodyfix.size(); i++)
				{
					vector<int> idx = { bodyfix[i] };
					vector<double> weight = { 1.0 };
					Eigen::Vector3d coord;
					coord = V6.row(bodyfix[i]).transpose();
					bodyfixidx.push_back(idx);
					bodyfixpos.push_back(coord);
					bodyfixweights.push_back(weight);
				}
				myH->p_DeformHandler->setconset(bodyfixidx, bodyfixweights, bodyfixpos);
				myH->p_DeformHandler->addconset(layer1_idx, layer1_weights, layer1_pos);
				//myH->p_DeformHandler->addconset(layer2_idx, layer2_weights, layer2_pos);
				//myH->p_DeformHandler->DeformInit1(v_id, V, V);
				myH->p_DeformHandler->econ_threshold = 2;
				start = clock();
				myH->p_DeformHandler->Deform(v_id, V7, V7);
				finish = clock();
				totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
				timeRecord << totaltime << " ";
				myH->p_DeformHandler->ModVMat(V7);

				Eigen::MatrixXd V8 = s * (fromsmpltoPc * V7.transpose());
				V8.colwise() += trans;
				Eigen::MatrixXd V9 = V8.transpose();
				sprintf_s(meshpath, 100, "data\\%d\\%d\\layer3.off", frame, i_th);
				fM->spCM->outMesh1(V9, meshpath);

				std::vector<int> bodyfix4 = { 3049, 4212, 1460, 4904, 4365, 877, 628, 1394, 5151, 4089, 7591, 1463, 1168, 4496, 4481, 1675, 1314, 4768, 5117 };
				std::vector<std::vector<int>> bodyfixidx4;
				std::vector<std::vector<double>> bodyfixweights4;
				std::vector<Eigen::Vector3d> bodyfixpos4;
				for (int i = 0; i < bodyfix4.size(); i++)
				{
					vector<int> idx = { bodyfix4[i] };
					vector<double> weight = { 1.0 };
					Eigen::Vector3d coord;
					coord = V7.row(bodyfix4[i]).transpose();
					bodyfixidx4.push_back(idx);
					bodyfixpos4.push_back(coord);
					bodyfixweights4.push_back(weight);
				}
				myH->p_DeformHandler->setconset(bodyfixidx, bodyfixweights, bodyfixpos);
				myH->p_DeformHandler->addconset(layer1_idx, layer1_weights, layer1_pos);
				myH->p_DeformHandler->addconset(layer2_idx, layer2_weights, layer2_pos);
				myH->p_DeformHandler->DeformInit1(v_id, V7, V7);
				myH->p_DeformHandler->Deform(v_id, V7, V7);
				myH->p_DeformHandler->ModVMat(V7);

				V8 = s * (fromsmpltoPc * V7.transpose());
				V8.colwise() += trans;
				V9 = V8.transpose();
				sprintf_s(meshpath, 100, "data\\%d\\%d\\layer4.off", frame, i_th);
				fM->spCM->outMesh1(V9, meshpath);
				/*if (newstate == 0)
				{
					sprintf_s(meshpath, 100, "data\\%d\\%d\\layer10.off", frame, i_th);
					fM->spCM->outMesh1(V9, meshpath);
				}*/


				//手部姿态
				PGMesh * pg1 = new PGMesh();
				outfile1 = ofstream("initialize.txt");
				for (int i = 0; i < 100; i++)
				{
					outfile1 << 0 << endl;
				}
				outfile1.close();
				sprintf_s(meshpath, 100, "data\\%d\\%d\\layer4.off", frame, i_th);
				OpenMesh::IO::read_mesh(*pg1, meshpath);
				Eigen::VectorXd tempV1 = fM->getVertexMatrix(pg1);
				delete pg1;
				Eigen::MatrixXd V10(3, fM->nV);
				for (int i = 0; i < fM->nV; i++)
				{
					V10.col(i) = tempV1.segment(3 * i, 3);
				}
				V10.colwise() -= trans;
				V10 = (1.0 / s) * V10;

				/*sRV2 = fromPctosmpl * sRV2;
				sRV = sRV2.transpose();*/
				//after
				Eigen::MatrixXd V11 = (fromPctosmpl * V10).transpose();
				//sRV = sRV2.transpose();

				fM->spCM->outMesh1(V11, "refMesh2.off");
				Eigen::MatrixXd V12(fM->nV, 3);
				myH = new mns::MyMeshHandler("refMesh2.off", "parameters3.txt", argc, argv);
				myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\male_hand\\DE\\newC.txt");
				myH->p_DeformHandler->DeformInit(v_id, V12, V12);
				std::vector<int> handfix = { 4087, 5177, 3038, 2896, 3523, 1678, 1668, 10450, 5221, 4840, 1709 };
				std::vector<std::vector<int>> handfixidx;
				std::vector<std::vector<double>> handfixweights;
				std::vector<Eigen::Vector3d> handfixpos;
				for (int i = 0; i < handfix.size(); i++)
				{
					vector<int> idx = { handfix[i] };
					vector<double> weight = { 1.0 };
					Eigen::Vector3d coord;
					coord = V11.row(handfix[i]).transpose();
					handfixidx.push_back(idx);
					handfixpos.push_back(coord);
					handfixweights.push_back(weight);
				}
				myH->p_DeformHandler->setconset(handfixidx, handfixweights, handfixpos);
				/*myH->p_DeformHandler->setconset(layer1_idx, layer1_weights, layer1_pos);
				myH->p_DeformHandler->addconset(layer2_idx, layer2_weights, layer2_pos);
				myH->p_DeformHandler->addconset(layer3_idx, layer3_weights, layer3_pos);*/
				myH->p_DeformHandler->addconset(fM->hand_feature_idx, fM->hand_feature_weights, fM->hand_feature_pos);

				//myH->p_DeformHandler->DeformInit1(v_id, V, V);
				start = clock();
				//if (newstate == 1)
				{
					myH->p_DeformHandler->Deform(v_id, V12, V12);
				}
				
				finish = clock();
				totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
				timeRecord << totaltime << endl;
				myH->p_DeformHandler->ModVMat(V12);

				/*Eigen::MatrixXd V13 = s * (fromsmpltoPc * V12.transpose());
				V13.colwise() += trans;
				Eigen::MatrixXd V14 = V13.transpose();*/
				Eigen::MatrixXd V14 = V12;
				sprintf_s(meshpath, 100, "data\\%d\\%d\\layer5.off", frame, i_th);
				fM->spCM->outMesh1(V14, meshpath);
				/*PGMesh *pg10 = new PGMesh();
				OpenMesh::IO::read_mesh(*pg10, meshpath);
				sprintf_s(meshpath, 100, "data\\%d\\%d\\layer5.off", frame, i_th);*/
			}
			
		}
		
	}

	//cout << y;
}

void testDenseCorresponding(int argc, char** argv)
{
	char *meshpath = new char[100];
	sprintf_s(meshpath, 100, "F:\\scan\\male_without_mouth\\0000.obj");
	PGMesh *pg = new PGMesh();
	OpenMesh::IO::read_mesh(*pg, meshpath);
	Eigen::VectorXd y(69);
	Eigen::VectorXd x(10);
	Eigen::VectorXd facecoef(200);
	Eigen::VectorXd gesturecoef(100);
	Eigen::VectorXd poseparm(400);
	Eigen::Matrix3d globalrotate, fromPctosmpl, fromsmpltoPc;
	double s;
	std::vector<int> layer = { 0, 1, 2, 3, 4, 5,6,7,8,9,10,11,12,13,14,15,16,17,18 };
	std::vector<int> transvec = { 0, 6, 7, 8, 12, 13, 14 };
	std::vector<int> transvec1 = { 5, 8, 11, 14 };
	std::vector<int> layer1 = { 3, 6, 9, 12 ,  4, 7, 10, 13, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28 };
	std::vector<int> layer2 = { 4, 7, 10, 13, 15, 16, 17, 18 };
	std::vector<int> layer3 = { 5, 8, 11, 14 };
	Eigen::Vector3d  trans, global_trans;
	FaceBodyModel * fM = new FaceBodyModel(10, 200, pg);
	
	Eigen::VectorXd out(fM->nV * 3);
	sprintf_s(meshpath, 100, "F:\\scan\\male_expressions\\DE\\newC.txt");
	fM->loadExpressionBases(meshpath, 200);
	sprintf_s(meshpath, 100, "F:\\project\\splocs-master\\newC.txt");
	fM->loadPoseBases(meshpath, 400);
	sprintf_s(meshpath, 100, "F:\\scan\\male_hand\\DE\\newC.txt");
	fM->loadgestureBases(meshpath, 100);
	sprintf_s(meshpath, 100, "F:\\scan\\smpl\\male_shape\\mean.txt");
	fM->loadMeanShape(meshpath);
	sprintf_s(meshpath, 100, "F:\\scan\\smpl\\male_shape\\shapebasis.txt");
	fM->loadShapeBases(meshpath);
	Eigen::MatrixXd deformedV(fM->nV, 3);
	fM->loadDenseCorrespondenceModule();
	std::vector<integer> v_id;
	PGMesh * tmpmesh = new PGMesh();
	fM->readBodyFeatureidx("body_feature_idx.txt");
	int a[19] = { 12, 15, 0, 16, 18, 20, 1, 4, 7, 17, 19, 21, 2, 5, 8, 24, 25, 26, 27 };
	fM->selectbodyfeature.resize(19);
	copy(a, a + 19, fM->selectbodyfeature.begin());
	fM->body_feature_idx1.resize(19);
	fM->body_feature_weights1.resize(19);
	for (int i = 0; i < 19; i++)
	{
		fM->body_feature_idx1[i] = fM->body_feature_idx[a[i]];
		fM->body_feature_weights1[i] = fM->body_feature_weights[a[i]];
	}
	for (int frame = 5361; frame < 5362; frame++)
	{
		fM->load_pointcloud(frame);
		for (int i_th = 0; i_th < 3; i_th++)
		{
			sprintf_s(meshpath, 100, "data\\%d\\%d\\coef.txt", frame, i_th);
			readcoefficient(meshpath, x, y);
			sprintf_s(meshpath, 100, "data\\%d\\%d\\global.txt", frame, i_th);
			read_global(meshpath, globalrotate, fromPctosmpl, fromsmpltoPc, s, trans, global_trans);
			sprintf_s(meshpath, 100, "data\\%d\\%d\\anchorpoint.txt", frame, i_th);
			//fM->spCM->loadAnchorData(meshpath);
			//fM->spCM->presolve();
			sprintf_s(meshpath, 100, "data\\%d\\%d\\bodyconset.txt", frame, i_th);
			fM->body_feature_pos1.resize(19);
			readBodyFeaturePos(meshpath, fM->body_feature_pos1);
			fM->body_feature_pos.resize(19);
			for (int i = 1; i < 19; i++)
			{
				fM->body_feature_pos[i] = fM->body_feature_pos1[i];
			}
			sprintf_s(meshpath, 100, "data\\%d\\%d\\handconset.txt", frame, i_th);
			fM->hand_feature_idx.resize(42);
			fM->hand_feature_weights.resize(42);
			fM->hand_feature_pos.resize(42);
			readHandFeaturePos(meshpath, fM->hand_feature_pos, fM->hand_feature_idx, fM->hand_feature_weights);
			sprintf_s(meshpath, 100, "data\\%d\\%d\\faceconset.txt", frame, i_th);
			fM->face_feature_idx.resize(54);
			fM->face_feature_pos.resize(54);
			fM->face_feature_weights.resize(54);
			readFaceFeaturePos(meshpath, fM->face_feature_pos, fM->face_feature_idx, fM->face_feature_weights);
			std::vector<std::vector<int>> denceidx;
			std::vector<std::vector<double>> denceweights;
			std::vector<Eigen::Vector3d> dencepos;

			std::vector<std::vector<int>> layer1_idx(layer1.size());
			std::vector<std::vector<double>> layer1_weights(layer1.size());
			std::vector<Eigen::Vector3d> layer1_pos(layer1.size());

			std::vector<std::vector<int>> layer2_idx(layer2.size());
			std::vector<std::vector<double>> layer2_weights(layer2.size());
			std::vector<Eigen::Vector3d> layer2_pos(layer2.size());

			std::vector<std::vector<int>> layer3_idx(layer3.size());
			std::vector<std::vector<double>> layer3_weights(layer3.size());
			std::vector<Eigen::Vector3d> layer3_pos(layer3.size());

			std::vector<std::vector<int>> facecontour_idx;
			std::vector<std::vector<double>> facecontour_weights;
			std::vector<Eigen::Vector3d> facecontour_pos;
			sprintf_s(meshpath, 100, "data\\%d\\%d\\facecontourconset.txt", frame, i_th);
			readFacecontourFeaturePos(meshpath, facecontour_pos, facecontour_idx, facecontour_weights);

			std::vector<int> layer4 = { 15, 16, 17, 18 };
			std::vector<std::vector<int>> layer4_idx(layer4.size());// +facecontour_pos.size());
			std::vector<std::vector<double>> layer4_weights(layer4.size());// +facecontour_pos.size());
			std::vector<Eigen::Vector3d> layer4_pos(layer4.size());// +facecontour_pos.size());

			for (int i = 0; i < layer1.size(); i++)
			{
				std::vector<int> idx = fM->body_feature_idx1[layer1[i]];
				std::vector<double> weight = fM->body_feature_weights1[layer1[i]];
				Eigen::Vector3d pos = fM->body_feature_pos[layer1[i]];
				layer1_idx[i] = idx;
				layer1_weights[i] = weight;
				layer1_pos[i] = pos;
			}
			for (int i = 0; i < layer2.size(); i++)
			{
				std::vector<int> idx = fM->body_feature_idx1[layer2[i]];
				std::vector<double> weight = fM->body_feature_weights1[layer2[i]];
				Eigen::Vector3d pos = fM->body_feature_pos[layer2[i]];
				layer2_idx[i] = idx;
				layer2_weights[i] = weight;
				layer2_pos[i] = pos;
			}
			for (int i = 0; i < layer3.size(); i++)
			{
				std::vector<int> idx = fM->body_feature_idx1[layer3[i]];
				std::vector<double> weight = fM->body_feature_weights1[layer3[i]];
				Eigen::Vector3d pos = fM->body_feature_pos[layer3[i]];
				layer3_idx[i] = idx;
				layer3_weights[i] = weight;
				layer3_pos[i] = pos;
			}

			for (int i = 0; i < layer4.size(); i++)
			{
				std::vector<int> idx = fM->body_feature_idx1[layer4[i]];
				std::vector<double> weight = fM->body_feature_weights1[layer4[i]];
				Eigen::Vector3d pos = fM->body_feature_pos[layer4[i]];
				layer4_idx[i] = idx;
				layer4_weights[i] = weight;
				layer4_pos[i] = pos;
			}

		
			//手部姿态
			PGMesh * pg1 = new PGMesh();
			ofstream outfile1;
			outfile1 = ofstream("initialize.txt");
			for (int i = 0; i < 400; i++)
			{
				outfile1 << 0 << endl;
			}
			outfile1.close();
			sprintf_s(meshpath, 100, "data\\%d\\%d\\layer5.off", frame, i_th);
			OpenMesh::IO::read_mesh(*pg1, meshpath);
			Eigen::VectorXd tempV1 = fM->getVertexMatrix(pg1);
			delete pg1;
			Eigen::MatrixXd V10(3, fM->nV);
			for (int i = 0; i < fM->nV; i++)
			{
				V10.col(i) = tempV1.segment(3 * i, 3);
			}
			V10.colwise() -= trans;
			V10 = (1.0 / s) * V10;

			/*sRV2 = fromPctosmpl * sRV2;
			sRV = sRV2.transpose();*/
			//after
			Eigen::MatrixXd V11 = (fromPctosmpl * V10).transpose();
			//sRV = sRV2.transpose();

			fM->spCM->outMesh1(V11, "refMesh2.off");
			Eigen::MatrixXd V12(fM->nV, 3);
			mns::MyMeshHandler *myH = new mns::MyMeshHandler("refMesh2.off", "parameters1.txt", argc, argv);
			myH->p_DeformHandler->ReceiveBasis2("C.txt");
			myH->p_DeformHandler->DeformInit(v_id, V12, V12);
			
			//myH->p_DeformHandler->setconset(handfixidx, handfixweights, handfixpos);
			myH->p_DeformHandler->setconset(layer1_idx, layer1_weights, layer1_pos);
			myH->p_DeformHandler->addconset(layer2_idx, layer2_weights, layer2_pos);
			myH->p_DeformHandler->addconset(layer3_idx, layer3_weights, layer3_pos);
			//myH->p_DeformHandler->addconset(fM->face_feature_idx, fM->face_feature_weights, fM->face_feature_pos);
			//myH->p_DeformHandler->addconset(fM->hand_feature_idx, fM->hand_feature_weights, fM->hand_feature_pos);

			//myH->p_DeformHandler->DeformInit1(v_id, V, V);
			myH->p_DeformHandler->econ_threshold = 1.0;
			//myH->p_DeformHandler->Deform(v_id, V12, V12);
			myH->p_DeformHandler->ModVMat(V12);
			myH->p_DeformHandler->addconset(layer4_idx, layer4_weights, layer4_pos);
			myH->p_DeformHandler->econ_threshold = 0.6;
			myH->p_DeformHandler->DeformInit1(v_id, V12, V12);
			//myH->p_DeformHandler->Deform(v_id, V12, V12);
			myH->p_DeformHandler->ModVMat(V12);
			delete myH;
			Eigen::MatrixXd V13 = s * (fromsmpltoPc * V12.transpose());
			V13.colwise() += trans;
			Eigen::MatrixXd V14 = V13.transpose();
			sprintf_s(meshpath, 100, "data\\%d\\%d\\layer6.off", frame, i_th);
			fM->spCM->outMesh1(V14, meshpath);

			//手部姿态
			pg1 = new PGMesh();
			outfile1 = ofstream("initialize.txt");
			for (int i = 0; i < 400; i++)
			{
				outfile1 << 0 << endl;
			}
			outfile1.close();
			sprintf_s(meshpath, 100, "data\\%d\\%d\\layer6.off", frame, i_th);
			OpenMesh::IO::read_mesh(*pg1, meshpath);
			tempV1 = fM->getVertexMatrix(pg1);
			delete pg1;
			for (int i = 0; i < fM->nV; i++)
			{
				V10.col(i) = tempV1.segment(3 * i, 3);
			}
			V10.colwise() -= trans;
			V10 = (1.0 / s) * V10;

			/*sRV2 = fromPctosmpl * sRV2;
			sRV = sRV2.transpose();*/
			//after
			V11 = (fromPctosmpl * V10).transpose();
			//sRV = sRV2.transpose();

			fM->spCM->outMesh1(V11, "refMesh2.off");
			myH = new mns::MyMeshHandler("refMesh2.off", "parameters2_1.txt", argc, argv);
			myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\male_expressions\\DE\\C.txt");
			myH->p_DeformHandler->DeformInit(v_id, V12, V12);
			std::vector<int> facefix = { 5298, 1864, 6448, 8387, 4782, 5105, 1651 };
			std::vector<std::vector<int>> facefixidx;
			std::vector<std::vector<double>> facefixweights;
			std::vector<Eigen::Vector3d> facefixpos;
			for (int i = 0; i < facefix.size(); i++)
			{
				vector<int> idx = { facefix[i] };
				vector<double> weight = { 1.0 };
				Eigen::Vector3d coord;
				coord = V11.row(facefix[i]).transpose();
				facefixidx.push_back(idx);
				facefixpos.push_back(coord);
				facefixweights.push_back(weight);
			}
			myH->p_DeformHandler->setconset(facefixidx, facefixweights, facefixpos);
			//myH->p_DeformHandler->setconset(handfixidx, handfixweights, handfixpos);
			//myH->p_DeformHandler->setconset(fM->body_feature_idx1, fM->body_feature_weights1, fM->body_feature_pos1);
			myH->p_DeformHandler->addconset(fM->face_feature_idx, fM->face_feature_weights, fM->face_feature_pos);
			//myH->p_DeformHandler->addconset(fM->hand_feature_idx, fM->hand_feature_weights, fM->hand_feature_pos);

			//myH->p_DeformHandler->DeformInit1(v_id, V, V);
			//myH->p_DeformHandler->Deform(v_id, V12, V12);
			myH->p_DeformHandler->ModVMat(V12);

			V13 = s * (fromsmpltoPc * V12.transpose());
			V13.colwise() += trans;
			V14 = V13.transpose();
			sprintf_s(meshpath, 100, "data\\%d\\%d\\layer7.off", frame, i_th);
			fM->spCM->outMesh1(V14, meshpath);

			//手部姿态
			pg1 = new PGMesh();


			sprintf_s(meshpath, 100, "data\\%d\\%d\\layer7.off", frame, i_th);
			OpenMesh::IO::read_mesh(*pg1, meshpath);
			tempV1 = fM->getVertexMatrix(pg1);
			delete pg1;
			for (int i = 0; i < fM->nV; i++)
			{
				V10.col(i) = tempV1.segment(3 * i, 3);
			}
			V10.colwise() -= trans;
			V10 = (1.0 / s) * V10;

			/*sRV2 = fromPctosmpl * sRV2;
			sRV = sRV2.transpose();*/
			//after
			V11 = (fromPctosmpl * V10).transpose();
			//sRV = sRV2.transpose();

			fM->spCM->outMesh1(V11, "refMesh2.off");
			myH = new mns::MyMeshHandler("refMesh2.off", "parameters3.txt", argc, argv);
			myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\male_hand\\DE\\newC.txt");
			myH->p_DeformHandler->DeformInit(v_id, V12, V12);

			//myH->p_DeformHandler->setconset(handfixidx, handfixweights, handfixpos);
			myH->p_DeformHandler->setconset(layer1_idx, layer1_weights, layer1_pos);
			myH->p_DeformHandler->addconset(layer2_idx, layer2_weights, layer2_pos);
			myH->p_DeformHandler->addconset(layer3_idx, layer3_weights, layer3_pos);
			myH->p_DeformHandler->addconset(fM->hand_feature_idx, fM->hand_feature_weights, fM->hand_feature_pos);

			//myH->p_DeformHandler->DeformInit1(v_id, V, V);
			//myH->p_DeformHandler->Deform(v_id, V12, V12);
			myH->p_DeformHandler->ModVMat(V12);

			V13 = s * (fromsmpltoPc * V12.transpose());
			V13.colwise() += trans;
			V14 = V13.transpose();
			sprintf_s(meshpath, 100, "data\\%d\\%d\\layer8.off", frame, i_th);
			fM->spCM->outMesh1(V14, meshpath);
			sprintf_s(meshpath, 100, "data\\%d\\%d\\layer9.off", frame, i_th);
			fM->spCM->outMesh1(V14, meshpath);


			for (int lll = 0; lll < 3; lll++)
			{
				fM->findCoressPairs(frame, i_th, 9);
				sprintf_s(meshpath, 100, "data\\%d\\%d\\coress.txt", frame, i_th);
				readCorespondence(meshpath, dencepos, denceidx, denceweights);
				for (int i = 0; i < dencepos.size(); i++)
				{
					dencepos[i] = (1.0 / s) * (dencepos[i] - trans);
					dencepos[i] = fromPctosmpl * dencepos[i];
				}
				pg1 = new PGMesh();
				sprintf_s(meshpath, 100, "data\\%d\\%d\\layer9.off", frame, i_th);
				OpenMesh::IO::read_mesh(*pg1, meshpath);
				Eigen::VectorXd tempV1 = fM->getVertexMatrix(pg1);
				delete pg1;
				Eigen::MatrixXd V10(3, fM->nV);
				for (int i = 0; i < fM->nV; i++)
				{
					V10.col(i) = tempV1.segment(3 * i, 3);
				}
				V10.colwise() -= trans;
				V10 = (1.0 / s) * V10;
				Eigen::MatrixXd V11 = (fromPctosmpl * V10).transpose();
				
				pg1 = new PGMesh();
				
				delete pg1;

				
				fM->spCM->outMesh1(V11, "refMesh2.off");
				sprintf_s(meshpath, 100, "refMesh2.off", frame, i_th);
				mns::MyMeshHandler *myH = new mns::MyMeshHandler(meshpath, "parameters4.txt", argc, argv);
				myH->p_DeformHandler->ReceiveBasis2("C.txt");
				
				myH->p_DeformHandler->DeformInit(v_id, V12, V12);

				//myH->p_DeformHandler->setconset(handfixidx, handfixweights, handfixpos);
				myH->p_DeformHandler->setconset(fM->body_feature_idx1, fM->body_feature_weights1, fM->body_feature_pos1);
				myH->p_DeformHandler->addconset(fM->face_feature_idx, fM->face_feature_weights, fM->face_feature_pos);
				myH->p_DeformHandler->addconset(fM->hand_feature_idx, fM->hand_feature_weights, fM->hand_feature_pos);
				myH->p_DeformHandler->addconset(denceidx, denceweights, dencepos);

				//myH->p_DeformHandler->DeformInit1(v_id, V, V);
				myH->p_DeformHandler->Deform(v_id, V12, V12);
				myH->p_DeformHandler->ModVMat(V12);


				Eigen::MatrixXd V13 = s * (fromsmpltoPc * V12.transpose());
				V13.colwise() += trans;
				Eigen::MatrixXd V14 = V13.transpose();
				sprintf_s(meshpath, 100, "data\\%d\\%d\\layer9.off", frame, i_th);
				fM->spCM->outMesh1(V14, meshpath);
				
			}

		}
	}
}

void testDenseShapeOpt(int argc, char ** argv)
{
	char *meshpath = new char[100];
	sprintf_s(meshpath, 100, "F:\\scan\\male_without_mouth\\0000.obj");
	PGMesh *pg = new PGMesh();
	OpenMesh::IO::read_mesh(*pg, meshpath);
	Eigen::VectorXd y(69);
	Eigen::VectorXd x(10);
	Eigen::VectorXd facecoef(200);
	Eigen::VectorXd gesturecoef(100);
	Eigen::VectorXd poseparm(400);
	Eigen::Matrix3d globalrotate, fromPctosmpl, fromsmpltoPc;
	double s;
	std::vector<int> layer = { 0, 1, 2, 3, 4, 5,6,7,8,9,10,11,12,13,14,15,16,17,18 };
	std::vector<int> transvec = { 0, 6, 7, 8, 12, 13, 14 };
	std::vector<int> transvec1 = { 5, 8, 11, 14 };
	std::vector<int> layer1 = { 3, 6, 9, 12 ,  4, 7, 10, 13, 15, 16, 17, 18, 5, 8, 11, 14 };
	std::vector<int> layer2 = { 4, 7, 10, 13, 15, 16, 17, 18 };
	std::vector<int> layer3 = { 5, 8, 11, 14 };
	Eigen::Vector3d  trans, global_trans;
	FaceBodyModel * fM = new FaceBodyModel(10, 200, pg);

	Eigen::VectorXd out(fM->nV * 3);
	sprintf_s(meshpath, 100, "F:\\scan\\male_expressions\\DE\\newC.txt");
	fM->loadExpressionBases(meshpath, 200);
	sprintf_s(meshpath, 100, "F:\\project\\splocs-master\\newC.txt");
	fM->loadPoseBases(meshpath, 400);
	sprintf_s(meshpath, 100, "F:\\scan\\male_hand\\DE\\newC.txt");
	fM->loadgestureBases(meshpath, 100);
	sprintf_s(meshpath, 100, "F:\\scan\\smpl\\male_shape\\mean.txt");
	fM->loadMeanShape(meshpath);
	sprintf_s(meshpath, 100, "F:\\scan\\smpl\\male_shape\\shapebasis.txt");
	fM->loadShapeBases(meshpath);
	Eigen::MatrixXd deformedV(fM->nV, 3);
	std::vector<integer> v_id;
	PGMesh * tmpmesh = new PGMesh();
	fM->loadDenseCorrespondenceModule();
	fM->readBodyFeatureidx("body_feature_idx.txt");
	int a[19] = { 12, 15, 0, 16, 18, 20, 1, 4, 7, 17, 19, 21, 2, 5, 8, 24, 25, 26, 27 };
	fM->selectbodyfeature.resize(19);
	copy(a, a + 19, fM->selectbodyfeature.begin());
	fM->body_feature_idx1.resize(19);
	fM->body_feature_weights1.resize(19);
	for (int i = 0; i < 19; i++)
	{
		fM->body_feature_idx1[i] = fM->body_feature_idx[a[i]];
		fM->body_feature_weights1[i] = fM->body_feature_weights[a[i]];
	}
	for (int frame = 5336; frame < 5337; frame++)
	{
		fM -> load_pointcloud(frame);
		for (int i_th = 0; i_th < 3; i_th++)
		{
			sprintf_s(meshpath, 100, "data\\%d\\%d\\coef.txt", frame, i_th);
			readcoefficient(meshpath, x, y);
			sprintf_s(meshpath, 100, "data\\%d\\%d\\global.txt", frame, i_th);
			read_global(meshpath, globalrotate, fromPctosmpl, fromsmpltoPc, s, trans, global_trans);
			sprintf_s(meshpath, 100, "data\\%d\\%d\\anchorpoint.txt", frame, i_th);
			//fM->spCM->loadAnchorData(meshpath);
			//fM->spCM->presolve();
			sprintf_s(meshpath, 100, "data\\%d\\%d\\bodyconset.txt", frame, i_th);
			fM->body_feature_pos1.resize(19);
			readBodyFeaturePos(meshpath, fM->body_feature_pos1);
			fM->body_feature_pos.resize(19);
			for (int i = 1; i < 19; i++)
			{
				fM->body_feature_pos[i] = fM->body_feature_pos1[i];
			}
			sprintf_s(meshpath, 100, "data\\%d\\%d\\handconset.txt", frame, i_th);
			fM->hand_feature_idx.resize(42);
			fM->hand_feature_weights.resize(42);
			fM->hand_feature_pos.resize(42);
			readHandFeaturePos(meshpath, fM->hand_feature_pos, fM->hand_feature_idx, fM->hand_feature_weights);
			sprintf_s(meshpath, 100, "data\\%d\\%d\\faceconset.txt", frame, i_th);
			fM->face_feature_idx.resize(54);
			fM->face_feature_pos.resize(54);
			fM->face_feature_weights.resize(54);
			readFaceFeaturePos(meshpath, fM->face_feature_pos, fM->face_feature_idx, fM->face_feature_weights);
			std::vector<std::vector<int>> denceidx;
			std::vector<std::vector<double>> denceweights;
			std::vector<Eigen::Vector3d> dencepos;
			fM->findCoressPairs(frame, i_th, 5);
			sprintf_s(meshpath, 100, "data\\%d\\%d\\coress.txt", frame, i_th);
			readCorespondence(meshpath, dencepos, denceidx, denceweights);
			for (int i = 0; i < dencepos.size(); i++)
			{
				dencepos[i] = (1.0 / s) * (dencepos[i] - trans);
				dencepos[i] = fromPctosmpl * dencepos[i];
			}
			//手部姿态
			PGMesh * pg1 = new PGMesh();
			ofstream outfile1;
			outfile1 = ofstream("initialize.txt");
			for (int i = 0; i < 400; i++)
			{
				outfile1 << 0 << endl;
			}
			outfile1.close();
			Eigen::MatrixXd V12(fM->nV, 3);
			
			for (int lll = 0; lll < 2; lll++)
			{
				pg1 = new PGMesh();
				sprintf_s(meshpath, 100, "data\\%d\\%d\\layer5.off", frame, i_th);
				OpenMesh::IO::read_mesh(*pg1, meshpath);
				Eigen::VectorXd tempV1 = fM->getVertexMatrix(pg1);
				delete pg1;
				Eigen::MatrixXd V10(3, fM->nV);
				for (int i = 0; i < fM->nV; i++)
				{
					V10.col(i) = tempV1.segment(3 * i, 3);
				}
				V10.colwise() -= trans;
				V10 = (1.0 / s) * V10;
				Eigen::MatrixXd V11 = (fromPctosmpl * V10).transpose();
				fM->spCM->setMesh(V11);
				fM->spCM->calculateFacesNormals();
				fM->spCM->calculateFacesFrame();
				fM->spCM->computeDiEdge();
				Eigen::VectorXd cm2 = fM->spCM->DiEdgeDataMatrix.row(0);
				pg1 = new PGMesh();
				sprintf_s(meshpath, 100, "data\\%d\\%d\\ref.off", frame, i_th);
				OpenMesh::IO::read_mesh(*pg1, meshpath);
				Eigen::VectorXd tempV2 = fM->getVertexMatrix(pg1);
				fM->spCM->setMesh(tempV2);
				fM->spCM->calculateFacesNormals();
				fM->spCM->calculateFacesFrame();
				fM->spCM->computeDiEdge();
				Eigen::VectorXd cm0 = fM->spCM->DiEdgeDataMatrix.row(0);
				delete pg1;

				Eigen::VectorXd cm1(2 * fM->nE);
				for (int i = 0; i < fM->nE; i++)
				{
					cm1(2 * i + 0) = cm2(2 * i + 0) - cm0(2 * i + 0);
					cm1(2 * i + 1) = cm2(2 * i + 1) / cm0(2 * i + 1) - 1;
				}
				/*sRV2 = fromPctosmpl * sRV2;
				sRV = sRV2.transpose();*/
				//after
				
				//sRV = sRV2.transpose();

				fM->spCM->outMesh1(V11, "refMesh2.off");
				sprintf_s(meshpath, 100, "data\\%d\\%d\\ref.off", frame, i_th);
				mns::MyMeshHandler *myH = new mns::MyMeshHandler(meshpath, "parameters6.txt", argc, argv);
				myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\smpl\\male_shape\\C.txt");
				myH->p_DeformHandler->cm1 = cm1;
				myH->p_DeformHandler->DeformInit5(v_id, V12, V12);

				//myH->p_DeformHandler->setconset(handfixidx, handfixweights, handfixpos);
				myH->p_DeformHandler->setconset(fM->body_feature_idx1, fM->body_feature_weights1, fM->body_feature_pos1);
				myH->p_DeformHandler->addconset(fM->face_feature_idx, fM->face_feature_weights, fM->face_feature_pos);
				myH->p_DeformHandler->addconset(fM->hand_feature_idx, fM->hand_feature_weights, fM->hand_feature_pos);
				myH->p_DeformHandler->addconset(denceidx, denceweights, dencepos);

				//myH->p_DeformHandler->DeformInit1(v_id, V, V);
				myH->p_DeformHandler->Deform1(v_id, V12, V12);
				myH->p_DeformHandler->ModVMat(V12);


				Eigen::MatrixXd V13 = s * (fromsmpltoPc * V12.transpose());
				V13.colwise() += trans;
				Eigen::MatrixXd V14 = V13.transpose();
				sprintf_s(meshpath, 100, "data\\%d\\%d\\layer5.off", frame, i_th);
				fM->spCM->outMesh1(V14, meshpath);
				fM->findCoressPairs(frame, i_th, 5);
				sprintf_s(meshpath, 100, "data\\%d\\%d\\coress.txt", frame, i_th);
				readCorespondence(meshpath, dencepos, denceidx, denceweights);
				for (int i = 0; i < dencepos.size(); i++)
				{
					dencepos[i] = (1.0 / s) * (dencepos[i] - trans);
					dencepos[i] = fromPctosmpl * dencepos[i];
				}
			}

			
		}
	}
}

wchar_t *GetWC(const char *c)
{
	const size_t cSize = strlen(c) + 1;
	size_t * ret = new size_t[cSize];
	wchar_t* wc = new wchar_t[cSize];
	mbstowcs_s(ret, wc, 100,  c, cSize);

	return wc;
}

void testpy()
{
	Py_SetPythonHome(GetWC("F:/Anaconda3/envs/python37"));
	Py_Initialize();//调用Py_Initialize()进行初始化
	PyRun_SimpleString("import sys");
	PyRun_SimpleString("sys.path.append('E:/Graphics/bachelor_thesis/Code/humanbody/humanbody/')");
	PyRun_SimpleString("sys.path.append('./')");
	PyObject * pModule = NULL;
	PyObject * pFunc = NULL;
	//PyRun_SimpleString("import tensorflow as tf");
	pModule = PyImport_ImportModule("TensorflowTest");//调用的Python文件名
	if (!pModule)
	{
		cout << "打开python文件失败";
		return;
	}
	pFunc = PyObject_GetAttrString(pModule, "Hello");//调用的函数名
	if (!pFunc)
	{
		cout << "无此方法";
		return;
	}

	////返回值
	PyObject *pReturn = NULL;
	pReturn = PyObject_CallObject(pFunc, NULL);//调用函数
	//										 //将返回值转换为int类型
	int result;
	PyArg_Parse(pReturn, "i", &result);//i表示转换成int型变量
	cout << "结果 = " << result << endl;

	////PyEval_CallObject(pFunc, NULL);//调用函数,NULL表示参数为空
	//Py_Finalize();//调用Py_Finalize,和Py_Initialize相对应的.
}



void getShapeBasis()
{
	char *meshpath = new char[100];
	sprintf_s(meshpath, 100, "F:\\scan\\smpl_H\\female_shape\\0000.obj");
	PGMesh *pg = new PGMesh();
	OpenMesh::IO::read_mesh(*pg, meshpath);
	AdaptiveSubdivision * adp = new AdaptiveSubdivision(pg);
	Eigen::MatrixXd ShapeBasis(10, adp->nV * 3);
	Eigen::VectorXd V0 = adp->getVertexVec(pg);

	for (int i = 0; i < 10; i++)
	{
		char *meshpath = new char[100];
		sprintf_s(meshpath, 100, "F:\\scan\\smpl_H\\female_shape\\%04d.obj", i+1);
		OpenMesh::IO::read_mesh(*pg, meshpath);
		ShapeBasis.row(i) = adp->getVertexVec(pg) - V0;
	}
	ofstream meantxt("F:\\scan\\smpl_H\\female_shape\\mean.txt");
	meantxt << V0;
	meantxt.close();
	ofstream basistxt("F:\\scan\\smpl_H\\female_shape\\shapebasis.txt");
	basistxt << ShapeBasis;
	basistxt.close();
}

void getAnchorPointsAndFaces()
{
	/*int ind[9] = { 4059, 598, 4786, 858, 1207, 3091, 8384, 5283, 1849 };
	int ind1[18] = { 331, 7591, 4096, 5094, 1873, 1677, 5099, 5318, 2867, 1702, 6823, 4505, 914, 1047, 999, 835, 4292, 4496 };
	int ind2[22] = { 6669, 6710, 3296, 3339, 3156, 6549, 5403, 5645, 5452, 2241, 2210, 2099, 2251, 2135, 2174, 2047, 5659, 5434, 5495, 5439, 6647, 3223 };
	*/
	int ind[12] = { 19522, 18827, 16651, 8026, 5785, 20689, 9815, 19705, 8834, 19903, 19270, 10080 };
	char *meshpath = new char[100];
	sprintf_s(meshpath, 100, "F:\\scan\\male_without_mouth\\0000.obj");
	PGMesh *pg = new PGMesh();
	OpenMesh::IO::read_mesh(*pg, meshpath);
	SparseLocalizedForConnectMap* ConnectMap = new SparseLocalizedForConnectMap(pg);
	ConnectMap->calculateFacesNormals();
	ConnectMap->calculateFacesFrame();
	ConnectMap->computeDiEdge();
	for (int i = 0; i < 4200; i++)
	{
		PGMesh *mesh2 = new PGMesh();
		char *path = new char[100];
		sprintf_s(path, 100, "F:\\scan\\male_without_mouth\\%04d.obj", i);
		OpenMesh::IO::read_mesh(*mesh2, path);
		char *path1 = new char[100];
		sprintf_s(path1, 100, "F:\\scan\\male_anchor_point\\%04d.txt", i);
		/*ofstream out(path1);
		for (int l = 0; l < 9; l++)
		{
			
			auto v = mesh2->vertex_handle(ind[l]);
			auto p = mesh2->point(v);
			out << p[0] << " " << p[1] << " " << p[2] << endl;
		}
		for (int l = 0; l < 18; l++)
		{
			
			auto v = mesh2->vertex_handle(ind1[l]);
			auto p = mesh2->point(v);
			out << p[0] << " " << p[1] << " " << p[2] << endl;
		}
		for (int l = 0; l < 22; l++)
		{
			
			auto v = mesh2->vertex_handle(ind2[l]);
			auto p = mesh2->point(v);
			out << p[0] << " " << p[1] << " " << p[2] << endl;
		}
		out.close();*/
		ConnectMap->updateDeformed1(mesh2);
		sprintf_s(path1, 100, "F:\\scan\\male_anchor_face\\%04d.txt", i);
		ofstream out1(path1);
		for (int j = 0; j < ConnectMap->faceAnchorFrames[1].size(); j++)
		{
			out1 << ConnectMap->faceAnchorFrames[1][j].localframe << endl;
		}
		out1.close();
	}
	
}

Eigen::MatrixXd getVertices(PGMesh *mesh_)
{
	Eigen::MatrixXd V(mesh_->n_vertices(), 3);
	for (auto vit = mesh_->vertices_begin(); vit != mesh_->vertices_end(); ++vit)
	{
		PGMesh::Point p = mesh_->point(*vit);
		
		V(vit->idx(), 0) = p[0];
		V(vit->idx(), 1)= p[1];
		V(vit->idx(), 2) = p[2];
	}
	return V;
}

void modifyVertices(Eigen::MatrixXd V, PGMesh *mesh_)
{
	for (auto vit = mesh_->vertices_begin(); vit != mesh_->vertices_end(); ++vit)
	{
		PGMesh::Point &p = mesh_->point(*vit);
		p[0] = V(vit->idx(), 0);
		p[1] = V(vit->idx(), 1);
		p[2] = V(vit->idx(), 2);
	}
}

void setSamplePoints1(string path)
{
	int ind[9] = { 4059, 598, 4786, 858, 1207, 3091, 8384, 5283, 1849 };
	int ind1[18] = { 331, 7591, 4096, 5094, 1873, 1677, 5099, 5318, 2867, 1702, 6823, 4505, 914, 1047, 999, 835, 4292, 4496 };
	int ind2[22] = { 6669, 6710, 3296, 3339, 3156, 6549, 5403, 5645, 5452, 2241, 2210, 2099, 2251, 2135, 2174, 2047, 5659, 5434, 5495, 5439, 6647, 3223 };
	ifstream input(path);
	std::vector<std::vector<double>> vertex_pos(9);
	std::vector<std::vector<double>> vertex_pos1(18);
	std::vector<std::vector<double>> vertex_pos2(22);
	for (int i = 0; i < 9; i++)
	{
		std::vector<double> a(3);
		input >> a[0];
		input >> a[1];
		input >> a[2];
		vertex_pos[i] = a;
	}
	for (int i = 0; i < 18; i++)
	{
		std::vector<double> a(3);
		input >> a[0];
		input >> a[1];
		input >> a[2];
		vertex_pos1[i] = a;
	}
	for (int i = 0; i < 22; i++)
	{
		std::vector<double> a(3);
		input >> a[0];
		input >> a[1];
		input >> a[2];
		vertex_pos2[i] = a;
	}
	input.close();
	ofstream out("conset.txt");
	for (int i = 0; i < 9; i++)
	{
		out << ind[i] << " ";
		
		out << vertex_pos[i][0] << " " << vertex_pos[i][1] << " " << vertex_pos[i][2] << endl;
	}
	for (int i = 2; i < 18; i++)
	{
		out << ind1[i] << " ";
		
		out << vertex_pos1[i][0] << " " << vertex_pos1[i][1] << " " << vertex_pos1[i][2] << endl;
	}
	for (int i = 0; i < 22; i++)
	{
		out << ind2[i] << " ";
		
		out << vertex_pos2[i][0] << " " << vertex_pos2[i][1] << " " << vertex_pos2[i][2] << endl;
	}
	out << -1;
	out.close();
	ofstream out0("conset0.txt");
	for (int i = 0; i < 9; i++)
	{
		out0 << ind[i] << " ";

		out0 << vertex_pos[i][0] << " " << vertex_pos[i][1] << " " << vertex_pos[i][2] << endl;
	}
	ofstream out1("conset1.txt");
	for (int i = 0; i < 18; i++)
	{
		out1 << ind1[i] << " ";

		out1 << vertex_pos1[i][0] << " " << vertex_pos1[i][1] << " " << vertex_pos1[i][2] << endl;
	}
	out1 << -1;
	out1.close();
	ofstream out2("conset2.txt");
	for (int i = 0; i < 22; i++)
	{
		out2 << ind2[i] << " ";

		out2 << vertex_pos2[i][0] << " " << vertex_pos2[i][1] << " " << vertex_pos2[i][2] << endl;
	}
	out2 << -1;
	out2.close();
}
void setSamplePoints(PGMesh *mesh_)
{
	/*int ind[13] = { 4867, 599, 4785, 5235, 1450, 2996, 3091, 2851, 6299, 5266, 1238, 1881, 3046 };
	int ind1[29] = {1655, 1663, 1910, 1732, 5175, 5063, 5064, 5185, 5136, 4851, 1429, 7449, 9555, 7476, 9772, 9176, 7590, 7396, 3041, 4479, 4471, 4465, 4458, 996, 1002, 1010, 1046, 1018,1016 };
	int ind2[18] = {5643, 5533, 5570,5633,2241,2000,2097,2136,6671,6555,6673,6738,3304,3406,3177,3299,3338,6686 };*/

	int ind[9] = { 4059, 598, 4786, 858, 1207, 3091, 8384, 5283, 1849 };
	int ind1[18] = { 331, 7591, 4096, 5094, 1873, 1677, 5099, 5318, 2867, 1702, 6823, 4505, 914, 1047, 999, 835, 4292, 4496 };
	int ind2[22] = { 6669, 6710, 3296, 3339, 3156, 6549, 5403, 5645, 5452, 2241, 2210, 2099, 2251, 2135, 2174, 2047, 5659, 5434, 5495, 5439, 6647, 3223 };
	int ind3[33] = { 5628, 5710, 5703, 5569, 5659, 5817, 5842, 5730, 5638, 6132,6150, 5913, 5973, 6070, 6025, 5597, 5705, 6031, 2287, 2133, 2270, 2173, 2193, 2617, 2636, 2499, 2516, 2367, 2408, 2296, 2168, 2717, 2288 };
	int ind4[42] = {9045, 8994, 1908, 6910, 6958, 7753, 8720, 9828, 7242, 6946, 8698, 9968, 3576, 5340, 6973, 6970, 9055, 8736, 6904, 331, 408, 7034, 436, 117, 7227, 7027, 7792, 9115, 9110, 9599, 9309, 9859, 9366, 9746, 9374, 9004, 8707, 57, 7052, 7281, 7465, 7644};
	ofstream out("conset.txt");
	for (int i = 0; i < 1; i++)
	{
		out <<ind4[ i] << " ";
		auto v = mesh_->vertex_handle(ind4[i]);
		auto p = mesh_->point(v);
		out << p[0] << " " << p[1] << " " << p[2] << endl;
	}
	/*for (int i = 0; i < 18; i++)
	{
		out << ind1[i] << " ";
		auto v = mesh_->vertex_handle(ind1[i]);
		auto p = mesh_->point(v);
		out<< p[0] << " " << p[1] << " " << p[2] << endl;
	}
	for (int i = 0; i < 22; i++)
	{
		out << ind2[i] << " ";
		auto v = mesh_->vertex_handle(ind2[i]);
		auto p = mesh_->point(v);
		out << p[0] << " " << p[1] << " " << p[2] << endl;
	}*/
	out << -1;
	out.close();
	ofstream out0("conset0.txt");
	for (int i = 0; i < 9; i++)
	{
		out0 << ind[i] << " ";
		auto v = mesh_->vertex_handle(ind[i]);
		auto p = mesh_->point(v);
		out0 << p[0] << " " << p[1] << " " << p[2] << endl;
	}
	ofstream out1("conset1.txt");
	for (int i = 0; i < 18; i++)
	{
		out1 << ind1[i] << " ";
		auto v = mesh_->vertex_handle(ind1[i]);
		auto p = mesh_->point(v);
		out1 << p[0] << " " << p[1] << " " << p[2] << endl;
	}
	out1 << -1;
	out1.close();
	ofstream out2("conset2.txt");
	for (int i = 0; i < 22; i++)
	{
		out2 << ind2[i] << " ";
		auto v = mesh_->vertex_handle(ind2[i]);
		auto p = mesh_->point(v);
		out2 << p[0] << " " << p[1] << " " << p[2] << endl;
	}
	out2 << -1;
	out2.close();
}

Eigen::VectorXd getDE(string path)
{
	ifstream input(path);
	Eigen::VectorXd x(65380);
	//Eigen::VectorXd x(41328);
	for (int i = 0; i < x.size(); i++)
	{
		input >> x(i);
	}
	input.close();
	return x;
}

void reconstruct_from_network()
{
	std::vector<int> valid_frame(3943);
	ifstream input("valid_frame.txt");
	for (int i = 0; i < valid_frame.size(); i++)
	{
		input >> valid_frame[i];
	}
	PGMesh *mesh2 = new PGMesh();
	OpenMesh::IO::read_mesh(*mesh2, "F:\\scan\\male_without_mouth\\0009.obj");
	//setSamplePoints1("anchor_point.txt");
	//setSamplePoints(mesh2);
	char *meshpath = new char[100];
	sprintf_s(meshpath, 100, "F:\\scan\\male_without_mouth\\0000.obj");
	PGMesh *pg = new PGMesh();
	OpenMesh::IO::read_mesh(*pg, meshpath);
	SparseLocalizedForConnectMap* ConnectMap = new SparseLocalizedForConnectMap(pg);
	ConnectMap->calculateFacesNormals();
	ConnectMap->calculateFacesFrame();
	ConnectMap->computeDiEdge();
	//ConnectMap->updateDeformed1(mesh2);
	int ref = 0;
	for (int frame = 0; frame < valid_frame.size(); frame++)
	{
		
		int currentframe = valid_frame[frame];
		int ref1 = (currentframe / 21) * 21;
		sprintf_s(meshpath, 100, "F:\\scan\\male_without_mouth\\%04d.obj", currentframe);
		OpenMesh::IO::read_mesh(*mesh2, meshpath);
		setSamplePoints(mesh2);
		if (ref1 != ref)
		{
			ref = ref1;
			sprintf_s(meshpath, 100, "F:\\scan\\male_without_mouth\\%04d.obj", ref1);
			OpenMesh::IO::read_mesh(*pg, meshpath);
			ConnectMap->setMesh(pg);
			ConnectMap->calculateFacesNormals();
			ConnectMap->calculateFacesFrame();
			ConnectMap->computeDiEdge();
		}
		sprintf_s(meshpath, 100, "face_anchor\\%04d.txt", frame);
		//ConnectMap->updateDeformed1(meshpath);
		ConnectMap->updateDeformed1(mesh2);
		ConnectMap->loadAnchorData("conset.txt");
		ConnectMap->presolve();
		Eigen::VectorXd l(2 * ConnectMap->nE);
		sprintf_s(meshpath, 100, "de\\%04d.txt", frame);
		Eigen::VectorXd l2 = getDE(meshpath);
		//Eigen::VectorXd l2 = ConnectMap->W_single.transpose() * ConnectMap->shapeBasis;
		Eigen::VectorXd l0 = ConnectMap->DiEdgeDataMatrix.row(0);
		for (int j = 0; j < ConnectMap->nE; j++)
		{
			l(2 * j + 0) = l0(2 * j + 0) + l2(2 * j + 0);
			l(2 * j + 1) = l0(2 * j + 1) * (1 + l2(2 * j + 1));
		}

		Eigen::VectorXd l1(3 * ConnectMap->nV);
		sprintf_s(meshpath, 100, "F:\\scan\\reconst_49\\%04d.off", currentframe);
		ConnectMap->reconstructionFromDiEdges(l, l1, 1, 0.5, 1);

		ConnectMap->outMesh1(l1, meshpath);
	}
}

void transferObj(int argc, char ** argv)
{
	char *meshpath = new char[100];
	sprintf_s(meshpath, 100, "F:\\scan\\faust\\tr_reg_000.ply");
	PGMesh *pg = new PGMesh();
	OpenMesh::IO::read_mesh(*pg, meshpath);
	SparseLocalizedForConnectMap* ConnectMap = new SparseLocalizedForConnectMap(pg);
	ConnectMap->calculateFacesNormals();
	ConnectMap->calculateFacesFrame();
	ConnectMap->computeDiEdge();
	ConnectMap->presolve();
	Eigen::VectorXd x0 = ConnectMap->DiEdgeDataMatrix.row(0);
	Eigen::VectorXd x = getDE("F:\\scan\\faust\\DE\\RLA.txt");
	Eigen::VectorXd l2(3 * ConnectMap->nV);

	ConnectMap->reconstructionFromDiEdges(x, l2, 1, 0.5, 1);
	ConnectMap->outMesh1(l2, "F:\\scan\\faust\\DE\\RLA.off");
}

void sythesisData(int argc, char ** argv)
{
	PGMesh *mesh2 = new PGMesh();
	OpenMesh::IO::read_mesh(*mesh2, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_B27\\layer2\\0806.off");
	char *meshpath = new char[100];
	sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_B27\\layer2\\0806.off");
	PGMesh *pg = new PGMesh();
	OpenMesh::IO::read_mesh(*pg, meshpath);
	SparseLocalizedForConnectMap* ConnectMap = new SparseLocalizedForConnectMap(pg);
	ConnectMap->calculateFacesNormals();
	ConnectMap->calculateFacesFrame();
	ConnectMap->computeDiEdge();
	ConnectMap->presolve();
	Eigen::VectorXd x0 = ConnectMap->DiEdgeDataMatrix.row(0);
	PGMesh *mesh3 = new PGMesh();
	OpenMesh::IO::read_mesh(*mesh3, "F:\\project\\humanbody\\humanbody\\data\\22782\\0\\initial_global.off");
	Eigen::VectorXd v0 = ConnectMap->getVertexMatrix(mesh3);
	ConnectMap->setMesh(v0);
	ConnectMap->calculateFacesNormals();
	ConnectMap->calculateFacesFrame();
	ConnectMap->computeDiEdge();
	Eigen::VectorXd x1 = ConnectMap->DiEdgeDataMatrix.row(0);
	OpenMesh::IO::read_mesh(*mesh3, "F:\\project\\humanbody\\humanbody\\data\\22782\\0\\layer1.off");
	Eigen::VectorXd v1 = ConnectMap->getVertexMatrix(mesh3);
	ConnectMap->setMesh(v1);
	ConnectMap->calculateFacesNormals();
	ConnectMap->calculateFacesFrame();
	ConnectMap->computeDiEdge();
	Eigen::VectorXd x2 = ConnectMap->DiEdgeDataMatrix.row(0);
	Eigen::VectorXd x = getRLA(x1, x2);
	Eigen::VectorXd out_x0 = getLAfromRLA(x0, x);
	Eigen::VectorXd l1(3 * ConnectMap->nV);

	ConnectMap->reconstructionFromDiEdges(out_x0, l1, 1, 0.5, 1);

	ConnectMap->outMesh1(l1, "out_0.off");

	OpenMesh::IO::read_mesh(*mesh3, "F:\\project\\humanbody\\humanbody\\data\\1618\\0\\layer4.off");
	Eigen::VectorXd v4 = ConnectMap->getVertexMatrix(mesh3);
	ConnectMap->setMesh(v4);
	ConnectMap->calculateFacesNormals();
	ConnectMap->calculateFacesFrame();
	ConnectMap->computeDiEdge();
	Eigen::VectorXd x4 = ConnectMap->DiEdgeDataMatrix.row(0);
	OpenMesh::IO::read_mesh(*mesh3, "F:\\project\\humanbody\\humanbody\\data\\1618\\0\\layer5.off");
	Eigen::VectorXd v5 = ConnectMap->getVertexMatrix(mesh3);
	ConnectMap->setMesh(v5);
	ConnectMap->calculateFacesNormals();
	ConnectMap->calculateFacesFrame();
	ConnectMap->computeDiEdge();
	Eigen::VectorXd x5 = ConnectMap->DiEdgeDataMatrix.row(0);
	Eigen::VectorXd x_new = getRLA(x4, x5);
	Eigen::VectorXd out_x1 = getLAfromRLA(out_x0, x_new);
	Eigen::VectorXd l2(3 * ConnectMap->nV);

	ConnectMap->reconstructionFromDiEdges(out_x1, l2, 1, 0.5, 1);

	ConnectMap->outMesh1(l2, "out_1.off");


}

void getInitialMesh()
{
	PGMesh *mesh2 = new PGMesh();
	OpenMesh::IO::read_mesh(*mesh2, "F:\\scan\\male_without_mouth\\0012.obj");
	//setSamplePoints1("anchor_point.txt");
	setSamplePoints(mesh2);
	char *meshpath = new char[100];
	sprintf_s(meshpath, 100, "F:\\scan\\male_without_mouth\\0000.obj");
	PGMesh *pg = new PGMesh();
	OpenMesh::IO::read_mesh(*pg, meshpath);
	SparseLocalizedForConnectMap* ConnectMap = new SparseLocalizedForConnectMap(pg);
	ConnectMap->calculateFacesNormals();
	ConnectMap->calculateFacesFrame();
	ConnectMap->computeDiEdge();
	ConnectMap->updateDeformed1(mesh2);
	//ConnectMap->updateDeformed1("face_anchor.txt");
	ConnectMap->loadAnchorData("conset.txt");
	
	ConnectMap->presolve();
	int n_comp = 400;
	ConnectMap->sparseLocalizedDecomp(n_comp);
	ConnectMap->getBasis("newC.txt", "W.txt");
	ConnectMap->setW_single("initialize.txt");
	//ConnectMap->setW_single("W_unknown.txt");
	Eigen::VectorXd l = ConnectMap->W_single.transpose() * ConnectMap->shapeBasis;
	Eigen::VectorXd l2 = getDE("de/0009t.txt");
	//cout << l2 << endl;
	//Eigen::VectorXd l2 = ConnectMap->W_single.transpose() * ConnectMap->shapeBasis;
	Eigen::VectorXd l0 = ConnectMap->DiEdgeDataMatrix.row(0);
	ofstream file2("true_cm0.txt");
	file2 << l0 << endl;
	for (int j = 0; j < ConnectMap->nE; j++)
	{
		l(2 * j + 0) =   l0(2 * j + 0) + l2(2 * j + 0);
		l(2 * j + 1) = l0(2 * j + 1) * (1 + l2(2 * j + 1)) ;
	}
	
	Eigen::VectorXd l1(3 * ConnectMap->nV);
	
	ConnectMap->reconstructionFromDiEdges(l, l1, 1, 0.5, 1);
	
	ConnectMap->outMesh1(l1, "initial_mesh_test.off");
}

void reconstructionFromNComp()
{
	std::vector<int> valid_frame(400);
	ifstream input("valid_frame.txt");
	for (int i = 0; i < valid_frame.size(); i++)
	{
		//input >> valid_frame[i];
		valid_frame[i] = i;
	}
	PGMesh *mesh2 = new PGMesh();
	OpenMesh::IO::read_mesh(*mesh2, "F:\\scan\\male_without_mouth\\0009.obj");
	//setSamplePoints1("anchor_point.txt");
	//setSamplePoints(mesh2);
	char *meshpath = new char[100];
	sprintf_s(meshpath, 100, "F:\\scan\\male_without_mouth\\0000.obj");
	PGMesh *pg = new PGMesh();
	OpenMesh::IO::read_mesh(*pg, meshpath);
	SparseLocalizedForConnectMap* ConnectMap = new SparseLocalizedForConnectMap(pg);
	ConnectMap->calculateFacesNormals();
	ConnectMap->calculateFacesFrame();
	ConnectMap->computeDiEdge();
	//ConnectMap->updateDeformed1(mesh2);
	//for (int i_th_comp = 0; i_th_comp < ncomps.size(); i_th_comp++)
	{
		int n_comp = 160;
		int nowcomp = ConnectMap->nSeq;
		ConnectMap->nSeq = 3000;
		ConnectMap->sparseLocalizedDecomp(n_comp);

		ConnectMap->getBasis("D:\\hand_reconstruct\\160\\C.txt", "D:\\hand_reconstruct\\160\\W.txt");
		ConnectMap->nSeq = nowcomp;
		int ref = 0;
		for (int frame = 0; frame < valid_frame.size(); frame++)
		{

			int currentframe = valid_frame[frame];
			int ref1 = (currentframe / 20) * 21;
			int ith = currentframe / 20;
			int jth = currentframe % 20 + 1;
			sprintf_s(meshpath, 100, "F:\\scan\\male_hand\\%04d.obj", 21 * ith + jth);
			OpenMesh::IO::read_mesh(*mesh2, meshpath);
			sprintf_s(meshpath, 100, "D:\\hand_reconstruct\\dataset\\mesh%04d.off", frame);
			//OpenMesh::IO::write_mesh(*mesh2, meshpath);
			setSamplePoints(mesh2);
			if (ref1 != ref)
			{
				ref = ref1;
				sprintf_s(meshpath, 100, "F:\\scan\\male_hand\\%04d.obj", ref1);
				OpenMesh::IO::read_mesh(*pg, meshpath);
				ConnectMap->setMesh(pg);
				ConnectMap->calculateFacesNormals();
				ConnectMap->calculateFacesFrame();
				ConnectMap->computeDiEdge();
			}
			//sprintf_s(meshpath, 100, "face_anchor\\%04d.txt", frame);
			//ConnectMap->updateDeformed1(meshpath);
			ConnectMap->updateDeformed1(mesh2);
			//ConnectMap->loadAnchorData("conset.txt");
			ConnectMap->presolve();
			Eigen::VectorXd l(2 * ConnectMap->nE);
			//sprintf_s(meshpath, 100, "E:\\reconstruction\\250\\work\\de\\%04d.txt", frame);
			//Eigen::VectorXd l2 = getDE(meshpath);
			Eigen::VectorXd l2 = ConnectMap->W.row(frame) * ConnectMap->shapeBasis;
			Eigen::VectorXd l0 = ConnectMap->DiEdgeDataMatrix.row(0);
			for (int j = 0; j < ConnectMap->nE; j++)
			{
				l(2 * j + 0) = l0(2 * j + 0) + l2(2 * j + 0);
				l(2 * j + 1) = l0(2 * j + 1) * (1 + l2(2 * j + 1));
			}

			Eigen::VectorXd l1(3 * ConnectMap->nV);
			sprintf_s(meshpath, 100, "D:\\hand_reconstruct\\160\\mesh%04d.off", frame);
			ConnectMap->reconstructionFromDiEdges(l, l1, 1, 0.5, 1);

			ConnectMap->outMesh1(l1, meshpath);
		}
	}

}

void editMesh(int argc, char** argv)
{
	
	mns::MyMeshHandler* myH = new mns::MyMeshHandler("F:\\scan\\male_without_mouth\\0000.obj", argc, argv);
	myH->p_DeformHandler->ReceiveBasis2();
	PGMesh *mesh2 = new PGMesh();
	OpenMesh::IO::read_mesh(*mesh2, "F:\\scan\\male_expressions\\0_22.obj");
	setSamplePoints(mesh2);
	
	myH->p_DeformHandler->current_str = "conset0.txt";
	mns::MyPointSet ps1 =  myH->p_DeformHandler->GetConSetByInput("conset0.txt");
	std::vector<integer> v_id;
	PGMesh *mesh1 = new PGMesh();
	OpenMesh::IO::read_mesh(*mesh1, "F:\\scan\\male_without_mouth\\0000.obj");
	OpenMesh::IO::write_mesh(*mesh1, "out_ini.off");
	Eigen::MatrixXd V0 (mesh1->n_vertices(), 3);
	myH->Deformation(v_id, V0, V0);
	/*mns::MyPointSet ps2 = myH->p_DeformHandler->GetConSetByInput("conset1.txt");
	myH->p_DeformHandler->current_str = "conset1.txt";
	myH->Deformation(v_id, V0, V0);*/
	/*myH->p_DeformHandler->current_str = "conset2.txt";
	myH->Deformation(v_id, V0, V0);*/
	modifyVertices(V0, mesh1);
	OpenMesh::IO::write_mesh(*mesh1, "outmesh.obj");
}

void outputEdgeNo(PGMesh * mesh_)
{
	std::vector<bool> edgeIn(mesh_->n_edges());
	for (int i = 0; i < mesh_->n_edges(); i++)
	{
		edgeIn[i] = false;
	}
	for (int i = 0; i < 1719; i++)
	{
		PGMesh::FaceHandle fit = mesh_->face_handle(faceSet[i]);
		//PGMesh::HalfedgeHandle heit;
		PGMesh::FaceHalfedgeIter fhit = mesh_->fh_begin(fit);
		//PGMesh::HalfedgeHandle h0 = fhit;
		PGMesh::EdgeHandle e0 = mesh_->edge_handle(*fhit);
		++fhit;
		//PGMesh::HalfedgeHandle h1 = fhit;
		PGMesh::EdgeHandle e1 = mesh_->edge_handle(*fhit);
		++fhit;
		//PGMesh::HalfedgeHandle h2 = fhit;
		PGMesh::EdgeHandle e2 = mesh_->edge_handle(*fhit);
		edgeIn[e0.idx()] = true;
		edgeIn[e1.idx()] = true;
		edgeIn[e2.idx()] = true;
	}
	ofstream out("edgeSet.txt");
	for (int i = 0; i < edgeIn.size(); i++)
	{
		if (edgeIn[i])
		{
			out << i << endl;
		}
	}
}
void processingDETrainingSets();

void show()
{
	char *meshpath = new char[100];
	sprintf_s(meshpath, 100, "F:\\scan\\male_without_mouth\\0000.obj");
	PGMesh *pg = new PGMesh();
	OpenMesh::IO::read_mesh(*pg, meshpath);
	SparseLocalizedForConnectMap* ConnectMap = new SparseLocalizedForConnectMap(pg);
	ConnectMap->calculateFacesNormals();
	ConnectMap->calculateFacesFrame();
	ConnectMap->computeDiEdge();
	ConnectMap->presolve();
	Eigen::VectorXd l0 = ConnectMap->DiEdgeDataMatrix.row(0);
	int n_comp = 500;
	ConnectMap->sparseLocalizedDecomp(n_comp);
	ConnectMap->getBasis("F:\\scan\\MFaust\\scripts\\male\\mixedDE\\newC.txt", "F:\\scan\\MFaust\\scripts\\male\\mixedDE\\W.txt");
	//ConnectMap->getBasis("F:\\scan\\male_expressions\\DE\\newC.txt", "W.txt");
	for (int i = 0; i < 500; i++)
	{
		
		char *compPath = new char[100];
		sprintf_s(compPath, 100, "F:\\scan\\MFaust\\scripts\\male\\mixedDE\\basis\\%04d.obj", i+1);
		PGMesh *pg1 = new PGMesh();
		OpenMesh::IO::read_mesh(*pg1, compPath);
		Eigen::VectorXd V(3 * pg1->n_vertices());
		
		Eigen::VectorXd l = 3 * ConnectMap->shapeBasis.row(i);
		Eigen::VectorXd l0 =  ConnectMap->DiEdgeDataMatrix.row(0);
		Eigen::VectorXd l9 = ConnectMap->shapeBasis.row(i);
		for (int j = 0; j < ConnectMap->nE; j++)
		{
			l(2 * j + 0) = l(2 * j + 0) + l0(2 * j + 0);
			//l(2 * j + 1) = l(2 * j + 1) + l0(2 * j + 1);
			l(2 * j + 1) = l0(2 * j + 1) * (1 +  l(2 * j + 1));
		}
		//l = l + l0;
		Eigen::VectorXd l1(3 * ConnectMap->nV);
		Eigen::VectorXd edgeError(ConnectMap->nE);
		for (int h = 0; h < ConnectMap->nE; h++)
		{
			edgeError(h) = l9.segment(2 * h, 2).norm();
		}
		ConnectMap->reconstructionFromDiEdges(l, l1, 0, 0.5, 0);
		ConnectMap->visualizeComponent(l1, edgeError, compPath);

	}
}
//void processingCMTrainingSets()
//{
//	typedef  OpenMesh::PolyMesh_ArrayKernelT<> PGMesh;
//	int totalPersons = 90;
//	int PoseNo = 20;
//#pragma omp parallel for
//	for (int i = 0; i < totalPersons; i++)
//	{
//		char *meshpath1 =new char[100];
//		sprintf_s(meshpath1, 100, "G:\\multipose\\female\\Object%04d.obj", 21 * i + 0);
//		PGMesh *pg0 = new PGMesh();
//		OpenMesh::IO::read_mesh(*pg0, meshpath1);
//		SparsePCAForShapeGrad* sparseLocal = new SparsePCAForShapeGrad(pg0);
//		for (int j = 0; j < PoseNo; j++)
//		{
//			PGMesh *pg1 = new PGMesh();
//			char *meshpath2 = new char[100];
//			sprintf_s(meshpath2, 100, "G:\\multipose\\female\\Object%04d.obj", 21 * i + j + 1);
//			
//			OpenMesh::IO::read_mesh(*pg1, meshpath2);
//			sparseLocal->setDeformedMesh(pg1);
//			sparseLocal->computeDG();
//			sparseLocal->computeCM();
//			char *path3 = new char[100];
//			sprintf_s(path3, 100, "G:\\multipose\\female\\CM\\%04d.txt", 21 * i + j + 1);
//			sparseLocal->writeCM(path3);
//		}
//	}
//}

//int main()
//{
//	//processingCMTrainingSets();
//	//typedef  OpenMesh::PolyMesh_ArrayKernelT<> PGMesh;
//	std::vector<PGMesh *> pg;
//	PGMesh *pg0 = new PGMesh();
//	PGMesh *pg1 = new PGMesh();
//	char *meshpath1= new char[100];
//	sprintf_s(meshpath1, 100, "D:\\bodydataset\\registered\\multi-pose-woman1\\Object%04d.obj", 21);
//	char *meshpath2 = new char[100];
//	sprintf_s(meshpath2, 100, "D:\\bodydataset\\registered\\multi-pose-woman1\\Object%04d.obj", 27);
//	OpenMesh::IO::read_mesh(*pg0, meshpath1);
//	OpenMesh::IO::read_mesh(*pg1, meshpath2);
//	//pg.push_back(pg0);
//	//outputEdgeNo(pg0);
//	//FaceBodyModel *fb_model = new FaceBodyModel(10, 200, pg0);
//	//fb_model->loadPoseBases("E:\\project\\splocs-master\\C.txt");
//	//fb_model->outputPoseBases();
//	SparsePCAForShapeGrad* sparseLocal = new SparsePCAForShapeGrad(pg0, pg1);
//	//sparseLocal->compute_V_inv();
//	//compute CM
//	sparseLocal->computeDG();
//	sparseLocal->computeCM();
//
//	//set the anchor point and anchorR
//	std::vector<int> fixPointIdx;
//	std::vector<Eigen::VectorXd> fixPointPos;
//	std::vector<int> fixFaceIdx;
//	std::vector<Eigen::MatrixXd> fixFaceR;
//	fixPointIdx.resize(5);
//	fixPointPos.resize(5);
//	fixFaceIdx.resize(50);
//	fixFaceR.resize(50);
//	for (int i = 0; i < 5; i++)
//	{
//		fixPointIdx[i] = i;
//		fixPointPos[i] = sparseLocal->SpatialTemData[i][1];
//	}
//	for (int i = 0; i < 50; i++)
//	{
//		fixFaceIdx[i] =   i;
//		fixFaceR[i] = sparseLocal->Rs[ i];
//	}
//	sparseLocal->setAnchorPoints(fixPointIdx, fixPointPos);
//	sparseLocal->setAnchorR(fixFaceIdx, fixFaceR);
//
//	sparseLocal->computeS_matrix();
//	sparseLocal->build_matrix();
//
//	std::vector<Eigen::MatrixXd> R_ij;
//	std::vector<Eigen::MatrixXd> S_i;
//	R_ij.resize(sparseLocal->nE);
//	S_i.resize(sparseLocal->nF);
//	sparseLocal->vecToMat(sparseLocal->VecCM, R_ij, S_i);
//	Eigen::VectorXd vert = sparseLocal->reconstructFromRS(sparseLocal->mesh_, R_ij, S_i);
//
//	sparseLocal->outMesh(vert, "reconst.off");
//
//	return 0;
//}
void IICsetK(int k, SparseLocalizedForConnectMap* ConnectMap)
{
	ConnectMap->sparseLocalizedDecomp(k);
	cout << k;

	char *tmp = new char[120];
	sprintf_s (tmp, 120, "E:\\project\\splocs-master\\C%d.txt", k);//C存放位置
	string str1 = tmp;
	sprintf_s (tmp, 120,  "E:\\project\\splocs-master\\W%d.txt", k);//W存放位置
	string str2 = tmp;
	ConnectMap->getBasis(str1, str2);
//#pragma omp parallel for
	for (int i = 3; i <=3; i++)
	{
		sprintf_s(tmp,120,  "mesh_%04d.off", 4);//重构存放位置
		str1 = tmp;
		ConnectMap->constructFrameFromComponents(i, k, str1);
	}
}

//int main(int argc, char** argv)
//{
//	QApplication a(argc, argv);
//	Main* main_window = new Main;
//	main_window->setAttribute(Qt::WA_DeleteOnClose);
//	main_window->exec();
//
//	//saveExpressionTensor();
//	//loadBlendshape("F:\\facewarehouse\\Tester_150\\Blendshape\\shape.bs");
//	//testFaceRegistritation(argc, argv);
//	//generateDatafromModel(argc, argv);
//	//testHumanModel8_1_woman(argc, argv);
//	//testHumanModel7(argc, argv);
//	//sythesisData(argc, argv);
//	//testHumanModel8_1(argc, argv);
//	//transferObj(argc, argv);
//	//reconstructionFromNComp();
//	//show();
//	//testDenseShapeOpt(argc, argv);
//	//generateTwistingPose(argc, argv);
//	//generateyuanqingtwistingpose(argc, argv);
//	//testReconstSeq(argc, argv);
//	//testDenseCorresponding(argc, argv);
//	//processingDETrainingSets1();
//	//testHumanModel1(argc, argv);
// 	//testpy();
//	//preSubdivisionDatasets();
//	//getShapeBasis();
//	//show();
//	//editMesh(argc, argv);
//	//processingDETrainingSets();
//	/*int totalMesh = 4200;
//	char *meshpath = new char[100];
//	sprintf_s(meshpath, 100, "F:\\scan\\male_without_mouth\\0000.obj");
//	PGMesh *pg = new PGMesh();
//	OpenMesh::IO::read_mesh(*pg, meshpath);
//	SparseLocalizedForConnectMap* ConnectMap = new SparseLocalizedForConnectMap(pg);
//	ConnectMap->subdivision();
//	ConnectMap->loadData();*/
//	//preSubdivisionDatasets();
//	//getInitialMesh();
//	//editMesh(argc, argv);
//	//reconstruct_from_network();
//	//getAnchorPointsAndFaces();
//	//processingDETrainingSets();
//	return 0;
//
//}

int main5()
{
	//ofstream l_stream("verts_samba.txt");
	std::vector<PGMesh *> pg;
	PGMesh *pg0 = new PGMesh();
	char *meshpath1= new char[100];
	sprintf_s(meshpath1, 100, "G:\\scan\\smoothed_female\\0_0_pose.obj");//原始序列第一帧
	//sprintf_s(meshpath1, 100, "D:\\bodydataset\\registered\\multi-pose-woman1\\Object%04d.obj", 21);
	OpenMesh::IO::read_mesh(*pg0, meshpath1);
	//for (int i = 0; i < 438; i++)
	//{
	//	PGMesh *pg1 = new PGMesh();
	//	sprintf_s(meshpath1, 100, "E:\\new_STED\\original\\cylinder_obj\\res_new_%03d.off", i);//原始序列第一帧
	//	char *meshpath2 = new char[100];
	//	OpenMesh::IO::read_mesh(*pg1, meshpath1);
	//	sprintf_s(meshpath2, 100, "E:\\new_STED\\original\\cylinder\\mesh_%04d.off", i);//原始序列第一帧
	//	OpenMesh::IO::write_mesh(*pg1, meshpath2);
	//}
	SparseLocalizedForConnectMap *ConnectMap = new SparseLocalizedForConnectMap(pg0);
	ConnectMap->calculateFacesNormals();
	ConnectMap->calculateFacesFrame();
	ConnectMap->computeDiEdge();
	//for (int i = 0; i < 175; i++)
	//{
	//	cout << i << endl;
	//	PGMesh *pg1 = new PGMesh();
	//	sprintf_s(meshpath1, 100, "E:\\new_STED\\original\\samba1\\mesh_%04d.off", i);//原始序列第一帧
	//	OpenMesh::IO::read_mesh(*pg1, meshpath1);
	//	ConnectMap->updateDeformed(pg1);
	//	ConnectMap->calculateFacesNormals();
	//	ConnectMap->calculateFacesFrame();
	//	ConnectMap->computeDiEdge();
	//	l_stream << ConnectMap->DiEdgeDataMatrix.row(1) << endl;
	//}
	//l_stream.close();
	ConnectMap->nSeq = 1356;
	ConnectMap->presolve();
	for (int i = 150; i <= 150; i += 10)
	{
		;
		IICsetK(i, ConnectMap);
	}

	
	
	return 0;
}

void preSubdivisionDatasets()
{
	char *meshpath = new char[100];
	sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_B24\\deformed\\0000.obj");
	PGMesh *pg = new PGMesh();
	OpenMesh::IO::read_mesh(*pg, meshpath);
	AdaptiveSubdivision * adp = new AdaptiveSubdivision(pg);
	char *headpath = new char[100];
	sprintf_s(headpath, 100, "head.txt");
	adp->getHeadVertices(headpath);

	adp->setSubdividedArea();
	adp->setSubdivisionCoeffs();
	adp->subdivision();
	int totalMesh = 484;
#pragma omp parallel for
	for (int i = 0; i < totalMesh; i++)
	{
		//if (!(i % 21))
		{
			//int i_th_person = i / 21;
			//int j_th_pose = i % 21;
			char *meshpath = new char[100];
			sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_B24\\deformed\\%04d.obj",  i);
			char *meshpath1 = new char[100];
			sprintf_s(meshpath1, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_B24\\init\\%04d.obj",i);
			PGMesh *pg1 = new PGMesh();
			OpenMesh::IO::read_mesh(*pg1, meshpath);
			Eigen::MatrixXd V0 = adp->getVertexMatrix(pg1);
			Eigen::MatrixXd V1 = adp->getSubdivisionPosition(V0);
			adp->outMesh(V1, meshpath1);
		}
	}

}

int main1()
{
	//show();
	/*int totalMesh = 4200;
	char *meshpath = new char[100];
	sprintf_s(meshpath, 100, "G:\\scan\\smoothed_female\\0_0_pose.obj");
	PGMesh *pg = new PGMesh();
	OpenMesh::IO::read_mesh(*pg, meshpath);
	SparseLocalizedForConnectMap* ConnectMap = new SparseLocalizedForConnectMap(pg);
	ConnectMap->subdivision();
	ConnectMap->loadData();*/
	//processingDETrainingSets();
	//preSubdivisionDatasets();
	/*char *meshpath = new char[100];
	sprintf_s(meshpath, 100, "G:\\generate_poses\\0_0_pose.obj");
	PGMesh *pg = new PGMesh();
	OpenMesh::IO::read_mesh(*pg, meshpath);
	AdaptiveSubdivision * adp = new AdaptiveSubdivision(pg);
	char *headpath = new char[100];
	sprintf_s(headpath, 100, "head.txt");
	adp->getHeadVertices(headpath);

	adp->setSubdividedArea();
	adp->setSubdivisionCoeffs();
	adp->subdivision();
	Eigen::MatrixXd V0 = adp->getVertexMatrix(pg);
	Eigen::MatrixXd V1 = adp->getSubdivisionPosition(V0);
	adp->outMesh(V1, "sub.obj");*/
	return 0; 
}

int main2()
{
	/*char *meshpath = new char[100];
	sprintf_s(meshpath, 100, "D:\\bodydataset\\registered\\multi-pose-woman1\\Object0000.obj");
	PGMesh *pg = new PGMesh();
	OpenMesh::IO::read_mesh(*pg, meshpath);
	outputEdgeNo(pg)*/;
	char *meshpath = new char[100];
	sprintf_s(meshpath, 100, "E:\\WYP\\multiposewoman_20180509\\0_0_pose.obj");
	PGMesh *pg = new PGMesh();
	OpenMesh::IO::read_mesh(*pg, meshpath);
	FaceBodyModel * fM = new FaceBodyModel(15, 100, pg);
	char *p1 = new char[100];
	sprintf_s(p1, 100, "E:\\WYP\\blendshapes.txt");
	fM->loadShapeBases(p1);
	p1 = new char[100];
	sprintf_s(p1, 100, "E:\\WYP\\mean_shape.txt");
	fM->loadMeanShape(p1);
	fM->outputBlendshapes();
	return 0;
}

void processingDETrainingSets1()
{
	//typedef  OpenMesh::PolyMesh_ArrayKernelT<> PGMesh;
	int totalMesh = 10;
	char *meshpath = new char[100];
	sprintf_s(meshpath, 100, "F:\\scan\\smpl\\male_shape\\0000.obj");
	PGMesh *pg = new PGMesh();
	OpenMesh::IO::read_mesh(*pg, meshpath);
	SparseLocalizedForConnectMap* ConnectMap = new SparseLocalizedForConnectMap(pg);
	ConnectMap->calculateFacesNormals();
	ConnectMap->calculateFacesFrame();
	ConnectMap->computeDiEdge();
	Eigen::VectorXd cm0= ConnectMap->DiEdgeDataMatrix.row(0);
	//#pragma omp parallel for
	sprintf_s(meshpath, 100, "F:\\scan\\smpl\\male_shape\\C.txt");
	ofstream op(meshpath, ios::out);
	int i_th, j_th;
	for (int i = 0; i < totalMesh; i++)
	{
		i_th = i + 1;
		char *meshpath1 = new char[100];
		
		sprintf_s(meshpath1, 100, "F:\\scan\\smpl\\male_shape\\%04d.obj", i_th);
		
		
		PGMesh *pg0 = new PGMesh();
		OpenMesh::IO::read_mesh(*pg0, meshpath1);
		ConnectMap->updateDeformed(pg0);
		ConnectMap->calculateFacesNormals();
		ConnectMap->calculateFacesFrame();
		ConnectMap->computeDiEdge();
		Eigen::VectorXd cm1 = ConnectMap->DiEdgeDataMatrix.row(1);
		op << cm1 - cm0 << endl;
	}
	op.close();
}

void processingDETrainingSets()
{
	//typedef  OpenMesh::PolyMesh_ArrayKernelT<> PGMesh;
	int totalMesh = 100;
	char *meshpath = new char[100];
	sprintf_s(meshpath, 100, "F:\\scan\\faust\\tr_reg_000.ply");
	PGMesh *pg = new PGMesh();
	OpenMesh::IO::read_mesh(*pg, meshpath);
	SparseLocalizedForConnectMap* ConnectMap = new SparseLocalizedForConnectMap(pg);
	//#pragma omp parallel for
	int i_th, j_th;
	for (int i = 0; i < totalMesh; i++)
	{
		i_th = i / 21;
		j_th = i % 21;
		char *meshpath1 = new char[100];
		//if (j_th == 0)
		{
			sprintf_s(meshpath1, 100, "F:\\scan\\faust\\tr_reg_%03d.ply", i);
		}
		//else
		{
		//	continue;
			//sprintf_s(meshpath1, 100, "F:\\scan\\male_hand\\%04d.obj", 21 * i_th+ j_th);
		}
		
		/*if (i > 3600 && i < 4001)
		{
			sprintf_s(meshpath1, 100, "G:\\multipose\\male\\add\\Object%04d.obj", i - 3601);
		}*/
		PGMesh *pg0 = new PGMesh();
		OpenMesh::IO::read_mesh(*pg0, meshpath1);
		ConnectMap->updateDeformed(pg0);
		ConnectMap->calculateFacesNormals();
		ConnectMap->calculateFacesFrame();
		ConnectMap->computeDiEdge();
		char *meshpath2 = new char[100];
		sprintf_s(meshpath2, 100, "F:\\scan\\faust\\DE\\%03d.txt", i);
		fstream op(meshpath2, ios::out);
		op << ConnectMap->DiEdgeDataMatrix.row(1) << endl;
		op.close();
	}
}

void testHumanModel7_(int argc, char** argv)
{
	char* meshpath = new char[100];
	sprintf_s(meshpath, 100, "F:\\scan\\male_without_mouth\\0000.obj");
	PGMesh* pg = new PGMesh();
	OpenMesh::IO::read_mesh(*pg, meshpath);
	Eigen::VectorXd y(69);
	Eigen::VectorXd x(10);
	Eigen::VectorXd facecoef(200);
	Eigen::VectorXd gesturecoef(100);
	Eigen::VectorXd poseparm(400);
	Eigen::Matrix3d globalrotate, fromPctosmpl, fromsmpltoPc;
	double s;
	std::vector<int> layer = { 0, 1, 2, 3, 4, 5,6,7,8,9,10,11,12,13,14,15,16,17,18,19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32 };
	std::vector<int> transvec = { 0, 6, 7, 8, 12, 13, 14 };
	std::vector<int> transvec1 = { 5, 8, 11, 14 };
	std::vector<int> layer1 = { 12, 13, 14, 16, 17, 6, 0, 1, 2, 4, 5, 18, 19, 24, 25, 26, 27, 28, 29, 30, 31, 32 };
	std::vector<int> layer2 = { 20, 21, 22, 23, 7, 8, 10, 11 };
	std::vector<int> layer3 = { 5, 8, 11, 14 };
	Eigen::Vector3d  trans, global_trans;
	FaceBodyModel* fM = new FaceBodyModel(10, 200, pg);
	//FaceBodyModel * fM = new FaceBodyModel(10, 250, pg);
	Eigen::VectorXd out(fM->nV * 3);
	sprintf_s(meshpath, 100, "F:\\scan\\male_expressions\\DE\\newC.txt");
	fM->loadExpressionBases(meshpath, 200);
	/*sprintf_s(meshpath, 100, "F:\\project\\splocs-master\\newC.txt");
	fM->loadPoseBases(meshpath, 400);*/
	/*sprintf_s(meshpath, 100, "F:\\scan\\MFaust\\scripts\\male\\mixedDE\\newC.txt");
	fM->loadPoseBases(meshpath, 500);*/
	sprintf_s(meshpath, 100, "F:\\scan\\male_without_mouth\\DE\\C.txt");//female->male
	fM->loadPoseBases(meshpath, 400);
	sprintf_s(meshpath, 100, "F:\\scan\\female_hand\\DE\\newC.txt");
	fM->loadgestureBases(meshpath, 100);
	sprintf_s(meshpath, 100, "F:\\scan\\smpl\\male_shape\\mean.txt");//female->male
	fM->loadMeanShape(meshpath);
	sprintf_s(meshpath, 100, "F:\\scan\\smpl\\male_shape\\shapebasis.txt");//female->male
	fM->loadShapeBases(meshpath);

	Eigen::MatrixXd deformedV(fM->nV, 3);
	fM->readBodyFeatureidx1("body_feature_idx1.txt");//赋值fM->body_feature_idx
	int a[33] = { 0, 1, 2, 3, 4,5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27, 28, 29, 30, 31, 32 };
	fM->selectbodyfeature.resize(33);
	copy(a, a + 33, fM->selectbodyfeature.begin());
	fM->body_feature_idx1.resize(33);
	fM->body_feature_weights1.resize(33);
	for (int i = 0; i < 33; i++)
	{
		fM->body_feature_idx1[i] = fM->body_feature_idx[a[i]];
		fM->body_feature_weights1[i] = fM->body_feature_weights[a[i]];
	}
	std::vector<int> bodyadd_conset_idx = { 570, 8359, 4152, 5207, 751, 2994, 4093, 677, 4097, 1861 };
	sprintf_s(meshpath, 100, "F:\\scan\\male_without_mouth\\0000.obj");//自己加的
	std::auto_ptr<mns::MyMeshHandler> p_meshHandler(new mns::MyMeshHandler(meshpath, fM->body_feature_idx1, fM->body_feature_weights1));
	Eigen::VectorXd tmpV(fM->nV * 3);
	std::vector<int> smplbody0 = { 3500, 4985, 912, 5333, 1294, 3022, 2967, 6427, 3107, 6530,1619, 1665, 5212, 5165, 4515, 1029, 4486, 1001 };
	std::vector<int> fbmmbody0 = { 4785, 4956, 912, 5306, 1294, 2995, 2940, 6373, 3080, 6476 ,1619, 1665, 5185, 5138, 4488, 1029, 4590, 1001 };
	std::vector<int> smplbody1 = { 3500, 4985, 912, 5333, 1294, 3022, 2967, 6427, 3107, 6530 };
	std::vector<int> fbmmbody1 = { 4785, 4956, 912, 5306, 1294, 2995, 2940, 6373, 3080, 6476 };
	std::vector<int> smplbody2 = { 1619, 1665, 5212, 5165, 4515, 1029, 4486, 1001 };
	std::vector<int> fbmmbody2 = { 1619, 1665, 5185, 5138, 4488, 1029, 4590, 1001 };
	std::vector<int> smplbody3 = { 1970, 1930, 5669, 5560, 6723, 6581, 3323, 3203 };
	std::vector<int> fbmmbody3 = { 1970, 1930, 5642, 5533, 6669, 6527, 3296, 3176 };
	PGMesh* tmpmesh = new PGMesh();
	vector<int> frameset;
	int startframe = 0;
	for (int i = 0; i < (375); i++)
	{
		frameset.push_back(i);
	}

	for (int ss = 0; ss < 1; ss++)//for (int ss = 0; ss < frameset.size(); ss++)
	{
		int frame = frameset[ss];
		vector<int> ith_set = { 0 };
		for (int sss = 0; sss < 1; sss++)//for (int sss = 0; sss < ith_set.size(); sss++)
		{
			for (int ii = 1; ii < 101; ii++) {
				int i_th = ith_set[sss];
				int refFrame;
				refFrame = frame;

				//F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\coef\\coef_test_initial_pc.txt
				sprintf_s(meshpath, 100, "E:\\Graphics\\bachelor_thesis\\Code\\smplify-x-master\\output_smpl_server_new\\results\\%02d_img\\coef.txt", ii);//%04d
				readcoefficient(meshpath, x, y);

				//sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\anchorpoint\\%04d.txt", frame);
				//fM->spCM->loadAnchorData(meshpath);
				//fM->spCM->presolve();
				sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\featurepoint\\%04d.txt", frame);
				fM->body_feature_pos1.resize(33);
				readBodyFeaturePos1(meshpath, fM->body_feature_pos1);
				fM->body_feature_pos.resize(33);
				for (int i = 1; i < 33; i++)
				{
					fM->body_feature_pos[i] = fM->body_feature_pos1[i];
				}
				sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\hand_conset\\%04d.txt", frame);
				readHandFeaturePos(meshpath, fM->hand_feature_pos, fM->hand_feature_idx, fM->hand_feature_weights);
				std::vector<std::vector<int>> feet_idx;
				std::vector<std::vector<double>> feet_weights;
				std::vector<Eigen::Vector3d> feet_pos;
				sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\feet_conset\\%04d.txt", frame);
				readHandFeaturePos(meshpath, feet_pos, feet_idx, feet_weights);
				/*sprintf_s(meshpath, 100, "data\\%d\\%d\\faceconset.txt", frame, i_th);
				fM->face_feature_idx.resize(54);
				fM->face_feature_pos.resize(54);
				fM->face_feature_weights.resize(54);
				readFaceFeaturePos(meshpath, fM->face_feature_pos, fM->face_feature_idx, fM->face_feature_weights);*/
				int b[54] = { 9160, 3030, 152, 7173, 9811, 9249, 7258, 8782, 444, 4250, 5242, 1871, 1574, 6445, 4778, 5237, 932, 4336, 4407, 1160, 3098, 1099, 4613, 4578, 3294, 1127, 1050, 870, 4313, 4510, 6551, 3131, 3146, 1216, 2952, 1909, 1976, 6303, 5141, 5401, 5499, 5619, 2794, 1569, 7629, 7690, 4759, 1244, 4167, 6256, 1435, 890, 3063, 4303 };


				std::vector<std::vector<int>> layer1_idx;
				std::vector<std::vector<double>> layer1_weights;
				std::vector<Eigen::Vector3d> layer1_pos;
				sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\layer1_fp\\%04d.txt", frame);
				readHandFeaturePos(meshpath, layer1_pos, layer1_idx, layer1_weights);

				std::vector<std::vector<int>> layer2_idx;//+ feet_idx.size());
				std::vector<std::vector<double>> layer2_weights;// +feet_idx.size());
				std::vector<Eigen::Vector3d> layer2_pos;// +feet_idx.size());
				sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\layer2_fp\\%04d.txt", frame);
				readHandFeaturePos(meshpath, layer2_pos, layer2_idx, layer2_weights);

				std::vector<std::vector<int>> layer3_idx;
				std::vector<std::vector<double>> layer3_weights;
				std::vector<Eigen::Vector3d> layer3_pos;
				sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\layer3_fp\\%04d.txt", frame);
				readHandFeaturePos(meshpath, layer3_pos, layer3_idx, layer3_weights);



				std::vector<int> layer4 = { 15, 16, 17, 18 };
				std::vector<std::vector<int>> layer4_idx(layer4.size());// +facecontour_pos.size());
				std::vector<std::vector<double>> layer4_weights(layer4.size());// +facecontour_pos.size());
				std::vector<Eigen::Vector3d> layer4_pos(layer4.size());

				poseparm = fM->getPoseParm(x, y);
				ofstream outfile1("initialize.txt");
				for (int i = 0; i < 250; i++)
				{
					outfile1 << 0 << endl;
				}
				for (int i = 0; i < 250; i++)
				{
					outfile1 << 0 << endl;
				}
				outfile1.close();
				//cout << poseparm << endl;
				Eigen::MatrixXd V = fM->f_x(x, poseparm);


				Eigen::VectorXd ref = fM->generateShapeWithoutPose(x);
				Eigen::MatrixXd refMatrix(fM->nV, 3);
				for (int i = 0; i < fM->nV; i++)
				{
					refMatrix.row(i) = ref.segment(3 * i, 3);
				}
				fM->spCM->outMesh1(ref, "refmesh.off");
				sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\%04d.ply", frame);//原bodyHands_REGISTRATIONS_B24改为bodyHands_REGISTRATIONS_A08
				OpenMesh::IO::read_mesh(*tmpmesh, meshpath);
				Eigen::MatrixXd smplposMatrix = fM->getVertexMatrix1(tmpmesh);


				Eigen::MatrixXd smplFeaturepos(smplbody0.size(), 3);
				for (int i = 0; i < smplbody0.size(); i++)//加入特征点的位置
				{
					smplFeaturepos.row(i) = smplposMatrix.row(smplbody0[i]);
				}
				Eigen::MatrixXd smplFeaturepos1(smplbody1.size(), 3);
				Eigen::MatrixXd smplFeaturepos2(smplbody2.size(), 3);
				Eigen::MatrixXd smplFeaturepos3(smplbody3.size(), 3);
				std::vector<std::vector<int>> fbmmbody1_idx(smplbody1.size());
				std::vector<std::vector<double>> fbmmbody1_weight(smplbody1.size());
				std::vector<std::vector<int>> fbmmbody2_idx(smplbody2.size());
				std::vector<std::vector<double>> fbmmbody2_weight(smplbody2.size());
				std::vector<std::vector<int>> fbmmbody3_idx(smplbody3.size());
				std::vector<std::vector<double>> fbmmbody3_weight(smplbody3.size());
				for (int i = 0; i < smplbody1.size(); i++)
				{
					smplFeaturepos1.row(i) = smplposMatrix.row(smplbody1[i]);
					std::vector<int> idx = { smplbody1[i] };
					std::vector<double> weight = { 1.0 };
					fbmmbody1_idx[i] = idx;
					fbmmbody1_weight[i] = weight;
				}
				for (int i = 0; i < layer2.size(); i++)
				{
					smplFeaturepos2.row(i) = smplposMatrix.row(smplbody2[i]);
					std::vector<int> idx = { smplbody2[i] };
					std::vector<double> weight = { 1.0 };
					fbmmbody2_idx[i] = idx;
					fbmmbody2_weight[i] = weight;
				}

				for (int i = 0; i < smplbody3.size(); i++)
				{
					smplFeaturepos3.row(i) = smplposMatrix.row(smplbody3[i]);
					std::vector<int> idx = { smplbody1[3] };
					std::vector<double> weight = { 1.0 };
					fbmmbody3_idx[i] = idx;
					fbmmbody3_weight[i] = weight;
				}
				//Eigen::VectorXd currentFeaturepos = fM->getBodyPosfromCurrentV(V, layer1);
				//Eigen::VectorXd currentFeaturepos = fM->getFacePosfromCurrentV(V);
				//优化全局旋转参数R,t
				//fM->findRigidTransform(currentFeaturepos, globalrotate, global_trans, smplFeaturepos);
				//fM->findRigidTransform(currentFeaturepos, globalrotate, global_trans, fM->face_feature_pos);
				Eigen::VectorXd currentFeaturepos = fM->getBodyPosfromCurrentV(V, layer);
				////			//优化全局旋转参数R,t
				fM->findRigidTransform(currentFeaturepos, globalrotate, global_trans, fM->body_feature_pos);
				refMatrix = (globalrotate * refMatrix.transpose()).transpose();
				refMatrix.rowwise() += global_trans.transpose();
				fM->spCM->outMesh1(refMatrix, "refmesh.off");
				fM->spCM->setMesh(refMatrix);
				fM->spCM->calculateFacesNormals();
				fM->spCM->calculateFacesFrame();
				fM->spCM->computeDiEdge();
				Eigen::VectorXd cm0 = fM->spCM->DiEdgeDataMatrix.row(0);


				for (int i = 0; i < layer4.size(); i++)
				{
					std::vector<int> idx = fM->body_feature_idx1[layer4[i]];
					std::vector<double> weight = fM->body_feature_weights1[layer4[i]];
					Eigen::Vector3d pos = fM->body_feature_pos[layer4[i]];
					layer4_idx[i] = idx;
					layer4_weights[i] = weight;
					layer4_pos[i] = pos;
				}


				Eigen::MatrixXd sRV = globalrotate * (V.transpose());
				sRV.colwise() += global_trans;
				Eigen::MatrixXd sRV1 = sRV.transpose();

				Eigen::MatrixXd sRV2 = sRV1.transpose();
				Eigen::MatrixXd sRV3 = sRV2.transpose();
				sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\init\\%04d.off", frame);
				fM->spCM->outMesh1(sRV3, meshpath);
				sprintf_s(meshpath, 100, "generatedata\\EHF\\%02d_mesh_after_groundtruth.off", ii);
				fM->spCM->outMesh1(sRV1, meshpath);//***输出****refMesh2.off
				p_meshHandler->fitmode = 0;


				mns::MyMeshHandler* myH = new mns::MyMeshHandler("refMesh2.off", "parameters3.txt", argc, argv);
				myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\male_hand\\DE\\newC.txt");
				std::vector<integer> v_id;



				std::vector<int> facefix = { 5298, 1864, 6448, 8387, 4782, 5105, 1651 };
				std::vector<std::vector<int>> facefixidx;
				std::vector<std::vector<double>> facefixweights;
				std::vector<Eigen::Vector3d> facefixpos;
				for (int i = 0; i < facefix.size(); i++)
				{
					vector<int> idx = { facefix[i] };
					vector<double> weight = { 1.0 };
					Eigen::Vector3d coord;
					coord = sRV1.row(facefix[i]).transpose();
					facefixidx.push_back(idx);
					facefixpos.push_back(coord);
					facefixweights.push_back(weight);
				}
				//第一层--脸部
				myH->p_DeformHandler->DeformInit(v_id, V, V);
				myH->p_DeformHandler->setconset(layer1_idx, layer1_weights, layer1_pos);
				myH->p_DeformHandler->addconset(layer2_idx, layer2_weights, layer2_pos);
				//myH->p_DeformHandler->addconset(layer3_idx, layer3_weights, layer3_pos);
				//myH->p_DeformHandler->DeformInit1(v_id, V, V);
				myH->p_DeformHandler->econ_threshold = 1.0;
				myH->p_DeformHandler->Deform(v_id, V, V);
				myH->p_DeformHandler->ModVMat(V);

				sRV1 = V;
				sRV2 = (sRV1.transpose());
				sRV3 = sRV2.transpose();
				sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\layer1\\%04d.off", frame);
				fM->spCM->outMesh1(sRV3, meshpath);

				myH->p_DeformHandler->DeformInit1(v_id, V, V);

				myH->p_DeformHandler->addconset(layer3_idx, layer3_weights, layer3_pos);
				myH->p_DeformHandler->econ_threshold = 1.0;
				myH->p_DeformHandler->Deform(v_id, V, V);
				myH->p_DeformHandler->ModVMat(V);
				sRV1 = V;
				sRV2 = (sRV1.transpose());
				sRV3 = sRV2.transpose();
				sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\layer2\\%04d.off", frame);
				fM->spCM->outMesh1(sRV3, meshpath);

				/* 第二层 身体*/
				PGMesh* pg0 = new PGMesh();
				//outfile1 = ofstream("initialize.txt");
				//for (int i = 0; i < 100; i++)
				//{
				//	outfile1 << 0 << endl;
				//}
				//outfile1.close();
				sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\layer2\\%04d.off", frame);
				OpenMesh::IO::read_mesh(*pg0, meshpath);
				Eigen::VectorXd tempV = fM->getVertexMatrix(pg0);
				delete pg0;
				Eigen::MatrixXd V5(3, fM->nV);
				for (int i = 0; i < fM->nV; i++)
				{
					V5.col(i) = tempV.segment(3 * i, 3);
				}

				/*sRV2 = fromPctosmpl * sRV2;
				sRV = sRV2.transpose();*/
				//after
				Eigen::MatrixXd V6 = (V5).transpose();
				//sRV = sRV2.transpose();

				fM->spCM->outMesh1(V6, "refMesh2.off");
				Eigen::MatrixXd V7(fM->nV, 3);
				myH = new mns::MyMeshHandler("refMesh2.off", "parameters3.txt", argc, argv);
				myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\male_hand\\DE\\newC.txt");
				myH->p_DeformHandler->DeformInit(v_id, V7, V7);
				std::vector<int> handfix = { 4087, 5177, 3038, 2896, 3523, 1678, 1668, 10450, 5221, 4840, 1709 };
				std::vector<std::vector<int>> handfixidx;
				std::vector<std::vector<double>> handfixweights;
				std::vector<Eigen::Vector3d> handfixpos;
				for (int i = 0; i < handfix.size(); i++)
				{
					vector<int> idx = { handfix[i] };
					vector<double> weight = { 1.0 };
					Eigen::Vector3d coord;
					coord = V6.row(handfix[i]).transpose();
					handfixidx.push_back(idx);
					handfixpos.push_back(coord);
					handfixweights.push_back(weight);
				}
				myH->p_DeformHandler->setconset(handfixidx, handfixweights, handfixpos);
				myH->p_DeformHandler->addconset(fM->hand_feature_idx, fM->hand_feature_weights, fM->hand_feature_pos);
				//myH->p_DeformHandler->addconset(layer2_idx, layer2_weights, layer2_pos);
				//myH->p_DeformHandler->DeformInit1(v_id, V, V);
				myH->p_DeformHandler->econ_threshold = 2;
				myH->p_DeformHandler->Deform(v_id, V7, V7);
				myH->p_DeformHandler->ModVMat(V7);

				Eigen::MatrixXd V8 = (V7.transpose());
				Eigen::MatrixXd V9 = V8.transpose();
				sprintf_s(meshpath, 100, "F:\\scan\\scan_smpl\\bodyHands_REGISTRATIONS_A08\\layer3\\%04d.off", frame);
				fM->spCM->outMesh1(V9, meshpath);

				//std::vector<int> bodyfix4 = { 3049, 4212, 1460, 4904, 4365, 877, 628, 1394, 5151, 4089, 7591, 1463, 1168, 4496, 4481, 1675, 1314, 4768, 5117 };
				//std::vector<std::vector<int>> bodyfixidx4;
				//std::vector<std::vector<double>> bodyfixweights4;
				//std::vector<Eigen::Vector3d> bodyfixpos4;
				//for (int i = 0; i < bodyfix4.size(); i++)
				//{
				//	vector<int> idx = { bodyfix4[i] };
				//	vector<double> weight = { 1.0 };
				//	Eigen::Vector3d coord;
				//	coord = V7.row(bodyfix4[i]).transpose();
				//	bodyfixidx4.push_back(idx);
				//	bodyfixpos4.push_back(coord);
				//	bodyfixweights4.push_back(weight);
				//}
				//myH->p_DeformHandler->setconset(bodyfixidx, bodyfixweights, bodyfixpos);
				//myH->p_DeformHandler->addconset(layer1_idx, layer1_weights, layer1_pos);
				//myH->p_DeformHandler->addconset(layer2_idx, layer2_weights, layer2_pos);
				//myH->p_DeformHandler->DeformInit1(v_id, V7, V7);
				//myH->p_DeformHandler->Deform(v_id, V7, V7);
				//myH->p_DeformHandler->ModVMat(V7);

				//V8 = s * (fromsmpltoPc * V7.transpose());
				//V8.colwise() += trans;
				//V9 = V8.transpose();
				//sprintf_s(meshpath, 100, "data\\%d\\%d\\layer4.off", frame, i_th);
				//fM->spCM->outMesh1(V9, meshpath);


				////手部姿态
				//PGMesh * pg1 = new PGMesh();
				//outfile1 = ofstream("initialize.txt");
				//for (int i = 0; i < 100; i++)
				//{
				//	outfile1 << 0 << endl;
				//}
				//outfile1.close();
				//sprintf_s(meshpath, 100, "data\\%d\\%d\\layer4.off", frame, i_th);
				//OpenMesh::IO::read_mesh(*pg1, meshpath);
				//Eigen::VectorXd tempV1 = fM->getVertexMatrix(pg1);
				//delete pg1;
				//Eigen::MatrixXd V10(3, fM->nV);
				//for (int i = 0; i < fM->nV; i++)
				//{
				//	V10.col(i) = tempV1.segment(3 * i, 3);
				//}
				//V10.colwise() -= trans;
				//V10 = (1.0 / s) * V10;

				///*sRV2 = fromPctosmpl * sRV2;
				//sRV = sRV2.transpose();*/
				////after
				//Eigen::MatrixXd V11 = (fromPctosmpl * V10).transpose();
				////sRV = sRV2.transpose();

				//fM->spCM->outMesh1(V11, "refMesh2.off");
				//Eigen::MatrixXd V12(fM->nV, 3);
				//myH = new mns::MyMeshHandler("refMesh2.off", "parameters3.txt", argc, argv);
				//myH->p_DeformHandler->ReceiveBasis2("F:\\scan\\male_hand\\DE\\newC.txt");
				//myH->p_DeformHandler->DeformInit(v_id, V12, V12);
				//std::vector<int> handfix = { 4087, 5177, 3038, 2896, 3523, 1678, 1668, 10450, 5221, 4840, 1709 };
				//std::vector<std::vector<int>> handfixidx;
				//std::vector<std::vector<double>> handfixweights;
				//std::vector<Eigen::Vector3d> handfixpos;
				//for (int i = 0; i < handfix.size(); i++)
				//{
				//	vector<int> idx = { handfix[i] };
				//	vector<double> weight = { 1.0 };
				//	Eigen::Vector3d coord;
				//	coord = V11.row(handfix[i]).transpose();
				//	handfixidx.push_back(idx);
				//	handfixpos.push_back(coord);
				//	handfixweights.push_back(weight);
				//}
				//myH->p_DeformHandler->setconset(handfixidx, handfixweights, handfixpos);
				///*myH->p_DeformHandler->setconset(layer1_idx, layer1_weights, layer1_pos);
				//myH->p_DeformHandler->addconset(layer2_idx, layer2_weights, layer2_pos);
				//myH->p_DeformHandler->addconset(layer3_idx, layer3_weights, layer3_pos);*/
				//myH->p_DeformHandler->addconset(fM->hand_feature_idx, fM->hand_feature_weights, fM->hand_feature_pos);

				////myH->p_DeformHandler->DeformInit1(v_id, V, V);
				//myH->p_DeformHandler->Deform(v_id, V12, V12);
				//myH->p_DeformHandler->ModVMat(V12);

				//Eigen::MatrixXd V13 = s * (fromsmpltoPc * V12.transpose());
				//V13.colwise() += trans;
				//Eigen::MatrixXd V14 = V13.transpose();
				//sprintf_s(meshpath, 100, "data\\%d\\%d\\layer5.off", frame, i_th);
				//fM->spCM->outMesh1(V14, meshpath);

			}
		}
	}
	cout << "finish" << endl;
}

#pragma execution_character_set("utf-8")
Main::Main(QWidget* parent) :
	QDialog(parent), ui(new Ui::Main)
{
	ui->setupUi(this);
}

Main::~Main()
{
	delete ui;
}

void Main::LoadPANOMAN()
{
	FILE* fp;
	fp = _popen("LoadPanoman.sh", "r");
	_pclose(fp);
	generateDatafromModel(this->argc, this->argv);
	testHumanModel7(this->argc, this->argv);
	//testconstraint(this->argc, this->argv);
	ui->label_panoman->setText("finish");
}

void Main::LoadOpenpose()
{
	FILE* fp;
	fp = _popen("LoadOpenpose.sh", "r");
	_pclose(fp);
	ui->label_openpose->setText("finish");
}

void Main::LoadSmplifyx()
{
	FILE* fp;
	fp = _popen("LoadSmplifyx.sh", "r");
	_pclose(fp);
	ui->label_smplifyx->setText("finish");
}

void Main::UploadPhoto()
{
	QString file_path = QFileDialog::getOpenFileName(this, "open file dialog", "C:/Users/Administrator/Desktop", "Txt files(*.jpg)");
	if (file_path.isEmpty())
		return;
	
	QStringList list = file_path.split("/");
	ui->label_photo->setText(list[list.size() - 1]);

	string str = string((const char*)file_path.toLocal8Bit());
	this->photo_path = str;
	ofstream os;
	os.open("./photo_path.txt");
	os << str;
	os.close();

	QImage* img = new QImage();
	img->load(file_path);
	QImage scaledimg;
	scaledimg = img->scaled(ui->label->width(), ui->label->height(), Qt::KeepAspectRatio);
	ui->label->setPixmap(QPixmap::fromImage(scaledimg));

}
