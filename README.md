# HMS JANUS

## Installation

Janus is written in go. You will need to make sure that go is installed (see you system details for how to do that).

Getting and building Janus is relatively straightforward. You can clone the repository with:

 ```bash
git clone git@github.com:xgrau/janus.git # or https://git.sr.ht/~hms/janus
```

This will make a `janus` directory. You can then go into the directory and:

```bash
mamba create -n janus conda-forge::go==1.16.5  conda-forge::nlopt==2.7.0
mamba activate janus

# navigate to the nlopt installation folder, give it execution permissions and finish nlopt install like so:
# nlopt: https://github.com/go-nlopt/nlopt
cd ../../go/pkg/mod/github.com/go-nlopt/nlopt
chmod -R 766 nlopt@v0.0.0-20230219125344-443d3362dcb5/
cd nlopt@v0.0.0-20230219125344-443d3362dcb5/
bash install_nlopt.sh
# You will need to enter your password (and if you have run this before, you may need to delete the nlopt*/build directory).
# But after you have completed this, rerun the build of Janus.

# navigate back to janus home folder
cd /home/xavi/Programes/janus
# explicitly set system compiler
export CC=/usr/bin/gcc
export CXX=/usr/bin/g++
# set nlopt flags
export CGO_CFLAGS="-I${CONDA_PREFIX}/include"
export CGO_LDFLAGS="-L${CONDA_PREFIX}/lib -lnlopt"

# build janus
go build -x -v janus.go

```

## Documentation

### Instructions

This guide assumes that you have already installed Janus using the information in the [README](../README.md). If you have done so then you should be able to get a help message that looks something like this if you run Janus without options:

```bash
↪ ./janus
  -aa string
    	amino acids (WAG, JTT, LG)
  -aicc
    	use AICc instead of BIC (BIC is default)
  -b	estimate branch lengths
  -cpuprofile string
    	write cpu profile to file
  -g	gamma rate variation
  -mc string
    	model conf (hmsh)
  -min int
    	minimum number of tips required (default 11)
  -rm
    	(for dna) estimate heterogeneous rate matrix along with the base comp
  -s string
    	seq filename
  -t string
    	tree filename
  -tc string
    	tree conf (hmsh)
  -ue
    	uncertainty existence
  -ul
    	uncertainty location
  -w int
    	number of threads (default 4)
```

If you aren't getting any options like this and instead you get an error or command not found, then there is something wrong with your compilation or installation.

Sometimes you will get an error after compiling and running an example, and this can be because of nlopt. To ensure that you have nlopt correctly installed, you can go to where it is installed in your system, for mine it is `godir/github.com/go-nlopt/nlopt`. Then you can run `bash install_nlopt.sh`. You will need to enter your password (and if you have run this before, you may need to delete the `nlopt*/build` directory). But after you have completed this, rerun the build of Janus.

### Basic usage

A typical usage of Janus would be something like this (using files from `example/1`, see below) `hmsj -s seqs -t seqs.treefile.rr`. This will default to running 4 threads (you can change with `-w #`) and not running any uncertainty analyses. A more complicated run might look like `hmsj -s seqs -t seqs.treefile.rr -ue -ul -w 2`. This is running, with 2 threads, the whole analysis along with two uncertainty analyses (location and existence).

### An example

For the first example, we will look at a single shift in base composition on a nucleotide dataset containing 100 tips. This is found in the `examples/1` directory. You can run the command `hmsj -s seqs -t seqs.treefile.rr`

This run takes, on my machine:

```bash
________________________________________________________
Executed in   93.50 secs    fish           external
   usr time  345.96 secs  676.00 micros  345.96 secs
   sys time    4.45 secs   99.00 micros    4.45 secs
```

The output looks something like this:

```bash
using BIC
lnL: -57716.29053200117
gtr: 57562.775343402325 [0.7190374005384313 0.6913800900893816 0.9482206735527946 0.8069245998621986 1.1529823745871557]
end bf: 57556.67068210288 [0.22295872415411455 0.24946560654489985 0.23727055800328986 0.29030511129769576]
lnL: -57556.67068210287
num of models: 19
starting IC: 116529.43119639708
final est
On 1/19 (116078.318822)
 --  115855.52801210285 -57209.35745703728 116529.43119639708
115855.52801210285 < 116529.43119639708 208
On 2/19 (116132.121993)
On 3/19 (116274.433182)
On 4/19 (116300.139624)
On 5/19 (116341.318492)
On 6/19 (116379.543475)
On 7/19 (116405.259556)
 --  115870.16425452748 -57206.31394533112 115855.52801210285
On 8/19 (116444.054514)
On 9/19 (116445.827246)
On 10/19 (116470.557774)
 --  115872.8067096821 -57207.635172908434 115855.52801210285
On 11/19 (116471.032944)
 --  115873.89988594392 -57208.18176103934 115855.52801210285
On 12/19 (116473.667552)
On 13/19 (116482.930205)
On 14/19 (116485.446872)
 --  115873.2543610251 -57207.85899857993 115855.52801210285
On 15/19 (116488.598598)
 --  115874.76566592157 -57208.61465102817 115855.52801210285
On 16/19 (116493.503301)
 --  115873.6941244232 -57208.07888027898 115855.52801210285
On 17/19 (116501.905157)
 --  115872.40325212215 -57207.43344412846 115855.52801210285
On 18/19 (116513.289550)
 --  115871.55057689504 -57207.0071065149 115855.52801210285
On 19/19 (116513.474312)
 --  115874.47694003869 -57208.470288086726 115855.52801210285
115855.52801210285 208
Final models
-------
0 [0.2573300113613611 0.2479903767548054 0.2534781918909178 0.24120141999291567]
1 [0.12548169168606033 0.24796690621960957 0.20850184418022866 0.4180495579141015]
(((((taxon_84:0.0291006494,taxon_85:0.0056149309):0.4036113799,taxon_86:0.4699520525):0.0243112286,(taxon_87:0.1084375121,(taxon_88:0.0183472754,taxon_89:0.02570298):0.0476988329):0.2577841562):0.4285237293,((taxon_90:0.4406361241,((taxon_91:0.1319465195,taxon_92:0.1344482614):0.0522201609,taxon_93:0.1568579937):0.1539140664):0.3315615,((taxon_94:0.2497188151,(taxon_95:0.0543828237,taxon_96:0.0340701483):0.0885822543):0.3719240183,((taxon_97:0.1631747274,taxon_98:0.115467999):0.1032014134,(taxon_99:0.2462715145,taxon_100:0.2158905238):0.0388554934):0.2547711584):0.1040213436):0.2181392776):0.0593798503,(((((taxon_47:0.0435440405,taxon_48:0.0115760878):0.5138725008,(((taxon_49:0.0942153608,taxon_50:0.1059623707):0.1501340392,(taxon_51:0.0446941024,taxon_52:0.0567720239):0.1983205307):0.0420993791,((taxon_53:0.0783998339,taxon_54:0.0967997073):0.20808815,(taxon_55:0.1831535007,(taxon_56:0.0796749637,taxon_57:0.058563633):0.1021705634):0.1079279307):0.0186338288):0.2023866902):0.049649391,(taxon_58:0.1973308949,(taxon_59:0.0080273346,taxon_60:0.0059718274):0.2702534226):0.3975385491):0.1077853753,(taxon_61:0.6754161061,(((taxon_62:0.1608252267,(taxon_63:0.1533677626,taxon_64:0.1349166138):0.0324405003):0.2041533384,(((taxon_65:0.0852483194,taxon_69:0.1005230339):0.0000020307,((taxon_66:0.0473797319,taxon_67:0.0425342714):0.0173658952,taxon_68:0.0797519741):0.0269253322):0.087341706,(taxon_70:0.0248201015,taxon_71:0.02260475):0.2043236238):0.1324358284):0.0264965906,(((taxon_72:0.047448998,taxon_73:0.0421624736):0.1175111587,((taxon_74:0.0208017452,(taxon_75:0.0046767083,taxon_76:0.0058193279):0.0349563162):0.021228853,taxon_77:0.0437590517):0.0811335747):0.0388168833,(((taxon_78:0.0629594581,(taxon_79:0.0828991742,taxon_80:0.0740049114):0.0194332432):0.075482515,(taxon_81:0.0980610457,taxon_82:0.099331305):0.0350954999):0.0287269833,taxon_83:0.1802378844):0.0051717463):0.2116192352):0.181207565):0.0565298722):0.2266683612,(taxon_46:0.5349338722,((((taxon_17:0.300287953,taxon_18:0.3051437913)1:0.03536188,((((taxon_19:0.1771757538,(taxon_20:0.0417172653,taxon_21:0.0469891702)1:0.1720212123)1:0.0055502434,(taxon_22:0.2290096366,taxon_23:0.1913317527)1:0.0100933191)1:0.012056141,((taxon_24:0.0546750884,((taxon_25:0.0102871829,taxon_26:0.0201830405)1:0.0163405041,taxon_27:0.0467005098)1:0.0000020282)1:0.1880740494,(taxon_28:0.0970519574,taxon_29:0.0811956978)1:0.1367139869)1:0.0165928314)1:0.0538715892,(((((taxon_30:0.0711609203,(taxon_31:0.0741974868,taxon_32:0.0635874436)1:0.0138807916)1:0.0516904969,((taxon_33:0.0063432992,taxon_34:0.0017052451)1:0.1068079359,taxon_35:0.0957524081)1:0.0300140345)1:0.009717226,(taxon_36:0.0936033877,((taxon_37:0.0169131338,taxon_38:0.0210030931)1:0.0305553407,(taxon_39:0.0263743579,taxon_40:0.0217939705)1:0.010650434)1:0.0651724385)1:0.0419165182)1:0.0159186035,taxon_41:0.1604140558)1:0.0467336324,(taxon_42:0.0135764785,taxon_43:0.0245202184)1:0.15188547)1:0.1231965527)1:0.0565830247)1:0.1137117741,(taxon_44:0.0112785915,taxon_45:0.039312505)1:0.4040823076)1:0.0473490049,(((taxon_9:0.2129934332,taxon_10:0.1724437235)1:0.0310794189,(taxon_11:0.1128102438,((taxon_12:0.0065593908,taxon_13:0.0084351081)1:0.0662196324,(taxon_14:0.0930739017,(taxon_15:0.0746619409,taxon_16:0.0555324563)1:0.0326609319)1:0.0134571383)1:0.0143544272)1:0.1077615404)1:0.0830495696,((((taxon_3:0.1146290849,taxon_4:0.1374712348)1:0.0440315864,(taxon_5:0.0368008802,taxon_6:0.0315520816)1:0.1384159395)1:0.0174304305,(taxon_7:0.0404071631,taxon_8:0.0374400606)1:0.1343242931)1:0.060603527,(taxon_1:0.2577438564,taxon_2:0.2445111358)1:0.0395700894)1:0.0401412591)1:0.1565276968)1:0.1014932017):0.4452885124):0.0593798503);
```

There is an additional outfile `seqs.treefile.rr.gophy.results.tre` that you can display in `Figtree`. If you open the file and color the branches by `model` you will see something like the figure below. You will see that the nested model is model *1* and the base is model *0*.

<img src="/~hms/janus/blob/master/examples/1/example.png" alt="example" width="300"/>

This is suggesting that there is a model shift at the clade that includes `taxon_1 and taxon_17`. And from our output we can see that this model has the following composition A:0.125 C:0.248 G:0.209 T:0.418 and the base model has the following composition A: 0.257 C: 0.248 G: 0.253 T: 0.241. The simulation that generated these data had the composition of 0.25 for all nucleotides for the base model and A:0.0814 C:0.234 G:0.184 T:0.501 for the nested model at the clade `taxon_1 and taxon_17`.

You will also notice that there are other options for display including `modelpar` which will give you the parameters for the model (these are also displayed at the end of the run).

We can run another analysis that includes uncertainty in the existence of shifts. To do this, we can run the analysis `hmsj -s seqs -t tree -ue`. We prefer these runs to those that look at uncertainty in the location of shifts because those are hard to interpret on larger and more complex trees.

The output will look relatively similar with the exception of the output tree, which will have information on the location of the shift and the existience of the shift. At the end of the run, you will see something like this

```bash
...
On 19/19 (116513.474310)
 --  115874.47697652284 -57208.4703063288 115855.52803318176
((((taxon_17,taxon_18),((((taxon_19,(taxon_20,taxon_21)),(taxon_22,taxon_23)),((taxon_24,((taxon_25,taxon_26),taxon_27)),(taxon_28,taxon_29))),(((((taxon_30,(taxon_31,taxon_32)),((taxon_33,taxon_34),taxon_35)),(taxon_36,((taxon_37,taxon_38),(taxon_39,taxon_40)))),taxon_41),(taxon_42,taxon_43)))),(taxon_44,taxon_45)),(((taxon_9,taxon_10),(taxon_11,((taxon_12,taxon_13),(taxon_14,(taxon_15,taxon_16))))),((((taxon_3,taxon_4),(taxon_5,taxon_6)),(taxon_7,taxon_8)),(taxon_1,taxon_2))))1 1
115855.52803318176 116574.4332537194 ( -57579.17171076403 ) 1
exists: 1 1
...
```

This is showing each shift in newick format and then the result from the "exists" test which is just looking at the IC of the shift being at that node vs not (`1` means it is strongly supported).

And if you open the `seqs.treefile.rr.gophy.results.tre` in `Figtree` again you will see some additional options.

### Flags

```bash
  -aa string
    	amino acids (WAG, JTT, LG)
```

If you have an amino acid dataset, you use this option and have these models that you can use.

```bash
  -aicc
    	use AICc instead of BIC (BIC is default)
```

If you want to use AICc instead of BIC for all IC calcs.

```bash
  -b	estimate branch lengths
```

_Very experimental._ You can try to reestimate branch lengths using the heterogeneous model.

```bash
  -cpuprofile string
    	write cpu profile to file
```

This is just for us to check performance things.

```bash
  -g	gamma rate variation
```

This is to allow for gamma rate variation in sites.

```bash
  -mc string
    	model conf (hmsh)
```

This is for us to test things between this and the hringhorni version.

```bash
  -min int
    	minimum number of tips required (default 11)
```

The minimum number of tips allowed for a shift.

```bash
  -rm
    	(for dna) estimate heterogeneous rate matrix along with the base comp
```

_Very experimental._ Estimate a heterogeneous rate matrix in addition to just the base comp. This is probably not advisable and won't be very accurate unless you have a lot of tips for each model shift.

```bash
  -s string
    	seq filename
```

The sequence data in fasta format.

```bash
  -t string
    	tree filename
```

The tree in newick format with branch lengths measured as substitutions per site.

```bash
  -tc string
    	tree conf (hmsh)
```

The tree output from hringhorni. This is for our testing.

```bash
  -ue
    	uncertainty existence
```

Include checking for the uncertainty of a shift existing. 1 means high and 0 low.

```bash
  -ul
    	uncertainty location
```

_Very experimental._ Include checking for the uncertainty of a shift location. 1 means high and 0 low.

```bash
  -w int
    	number of threads (default 4)
```

Janus is parallel, so how many threads? Probably not more than you have on your machine.

### A note about optimization

Optimization of many parameters in a complex phylogenetic space is difficult. Janus is no different. We have attempted to develop the program so that it doesn't get stuck in local optima and so that it doesn't suffer from issues where it doesn't find the ML. The problem is exacerbated by likelihood space that is relatively rough and uneven. Think about using the hringhorni version if you are having particular problems.

### A note about hringhorni edition

In part to address the issues noted above and also just in the interest of speed (important for larger analyses), we have developed a version of this procedure and program written in `C` that is part of the hringhorni package. You can find that [here](http://git.sr.ht/~hms/hringhorni). This version of janus was the first and continues to be maintained as we try out new procedures here first.
