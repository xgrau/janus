package main

import (
	"bufio"
	"bytes"
	"flag"
	"fmt"
	"io"
	"log"
	"math"
	"os"
	"runtime/pprof"
	"sort"
	"strconv"
	"strings"
	"time"

	"gonum.org/v1/gonum/mat"

	"github.com/FePhyFoFum/gophy"
	"gonum.org/v1/gonum/floats"

	"golang.org/x/exp/rand"
)

// this will determine how many tips are not already covered by another model
func getVisible(nd *gophy.Node) (numvisible int) {
	dontinclude := map[string]bool{}
	for _, i := range nd.PreorderArray() {
		if _, ok := i.IData["shift"]; ok {
			for _, j := range i.GetTips() {
				dontinclude[j.Nam] = true
			}
		}
	}
	numvisible = len(nd.GetTips()) - len(dontinclude)
	return
}

type keyVal struct {
	Key   *gophy.Node
	Value float64
}

// returns a list with the smallest values
func sortAicMap(nodevalues map[*gophy.Node]float64) (ss []keyVal) {
	//sort the values
	ss = []keyVal{}
	for k, v := range nodevalues {
		ss = append(ss, keyVal{k, v})
	}

	sort.Slice(ss, func(i, j int) bool {
		return ss[i].Value < ss[j].Value
	})
	return
}

// mark the model
func markModel(nd *gophy.Node, modelnum int) {
	nd.IData["shift"] = modelnum
	nd.Nam = strconv.Itoa(modelnum)
	for _, i := range nd.PreorderArray() {
		if _, ok := i.IData["shift"]; !ok {
			i.IData["shift"] = modelnum
		}
	}
}

// mark the model
func unmarkModel(nd *gophy.Node, modelnum int) {
	delete(nd.IData, "shift")
	nd.Nam = ""
	for _, i := range nd.PreorderArray() {
		if _, ok := i.IData["shift"]; ok {
			if i.IData["shift"] == modelnum {
				delete(i.IData, "shift")
			}
		}
	}
}

func moveMarksToLabels(nd *gophy.Node) {
	for _, i := range nd.PreorderArray() {
		if len(i.Chs) == 0 {
			continue
		}
		if _, ok := i.IData["shift"]; ok {
			i.Nam = strconv.Itoa(i.IData["shift"])
		}
	}
}

func getNodeModels(curmodels []*gophy.DiscreteModel, curnodemodels map[*gophy.Node]int,
	y *gophy.DiscreteModel, nd *gophy.Node) (models []*gophy.DiscreteModel, nodemodels map[*gophy.Node]int, modelint int) {
	modelint = len(curmodels)
	markModel(nd, modelint)
	models = []*gophy.DiscreteModel{}
	for _, i := range curmodels {
		models = append(models, i)
	}
	models = append(models, y)
	nodemodels = make(map[*gophy.Node]int)
	for i, j := range curnodemodels {
		nodemodels[i] = j
	}
	for _, i := range nd.PreorderArray() {
		if i.IData["shift"] == modelint {
			nodemodels[i] = modelint
		}
	}
	return
}

// this is meant to check, after a deep node has been supported with a model shift, whether this is a nested one back to the revert
// must meet minimum requirements
func checkNested() {

}

func main() {
	rand.Seed(uint64(time.Now().UTC().UnixNano()))
	tfn := flag.String("t", "", "tree filename")
	afn := flag.String("s", "", "seq filename")
	uc1 := flag.Bool("ue", false, "uncertainty existence")
	uc2 := flag.Bool("ul", false, "uncertainty location")
	aa := flag.String("aa", "", "amino acids (WAG, JTT, LG)")
	bl := flag.Bool("b", false, "estimate branch lengths")
	gm := flag.Bool("g", false, "gamma rate variation")
	rm := flag.Bool("rm", false, "(for dna) estimate heterogeneous rate matrix along with the base comp")
	mintest := flag.Int("min", 11, "minimum number of tips required") //
	aicc := flag.Bool("aicc", false, "use AICc instead of BIC (BIC is default)")
	wks := flag.Int("w", 4, "number of threads")
	readsimpconf1 := flag.String("tc", "", "tree conf (hmsh)")
	readsimpconf2 := flag.String("mc", "", "model conf (hmsh)")
	cpuprofile := flag.String("cpuprofile", "", "write cpu profile to file")
	flag.Parse()
	minset := *mintest - 1 // the test, tests for > minset
	if len(*tfn) == 0 {
		flag.PrintDefaults()
		os.Exit(1)
	}
	if len(*afn) == 0 {
		flag.PrintDefaults()
		os.Exit(1)
	}
	if *cpuprofile != "" {
		f, err := os.Create(*cpuprofile)
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

	//
	var icfun func(float64, float64, int) (x float64)
	icfun = gophy.CalcBIC
	if *aicc {
		fmt.Fprintln(os.Stderr, "using AICc")
		icfun = gophy.CalcAICC
	} else {
		fmt.Fprintln(os.Stderr, "using BIC")
	}
	amino := false
	if len(*aa) > 0 {
		amino = true
		if *aa != "WAG" && *aa != "JTT" && *aa != "LG" {
			fmt.Fprintln(os.Stderr, "amino acid model must be WAG, JTT, or LG")
			os.Exit(1)
		}
	}

	if *bl {
		fmt.Fprintln(os.Stderr, "reestimating branch lengths with heterogeneous models")
	}

	if *gm {
		fmt.Fprintln(os.Stderr, "estimating with gamma rate variation")
	}

	//read tree
	t := gophy.ReadTreeFromFile(*tfn)
	seqs, patternsint, nsites, bf := gophy.ReadPatternsSeqsFromFile(*afn, !(amino))
	var patternval []float64
	if amino {
		patternval, _ = gophy.PreparePatternVecsProt(t, patternsint, seqs)
	} else {
		patternval, _ = gophy.PreparePatternVecs(t, patternsint, seqs)
	}
	modelparams := make([]float64, 5)
	for i := range modelparams {
		modelparams[i] = 1.0
	}
	//root model
	var M gophy.DiscreteModel

	if amino {
		x := gophy.NewProteinModel()
		x.M.SetBaseFreqs(bf)
		if *aa == "JTT" {
			x.SetRateMatrixJTT()
			x.M.SetBaseFreqs(bf)
			fmt.Println(x.M.BF)
		} else if *aa == "LG" {
			x.SetRateMatrixLG()
			x.M.SetBaseFreqs(bf)
			fmt.Println(x.M.BF)
		} else if *aa == "WAG" {
			x.SetRateMatrixWAG()
			x.M.SetBaseFreqs(x.M.BF)
			fmt.Println(x.M.BF)
		}
		x.M.SetupQGTR()
		M = x.M
	} else {
		x := gophy.NewDNAModel()
		x.M.SetBaseFreqs(bf)
		x.M.SetRateMatrix(modelparams) // equal rates
		x.M.SetupQGTR()
		M = x.M
	}
	if *gm {
		M.GammaNCats = 4
		M.GammaAlpha = 1.0
		M.GammaCats = gophy.GetGammaCats(M.GammaAlpha, M.GammaNCats, false)
	}

	l, useLog := initOptimization(t, &M, patternval, *wks, amino, &modelparams, *bl, *gm)

	//read an existing configuration file / probably for testing uncertainty

	var curmodels []*gophy.DiscreteModel
	var curnodemodels map[*gophy.Node]int
	var modelnodes []*gophy.Node
	var curparams float64
	var currentaic float64

	if len(*readsimpconf1) != 0 {
		if len(*readsimpconf2) == 0 {
			fmt.Println("Optimising parameters and calculating likelihood of prespecified configuration")
			t = gophy.ReadTreeFromFile(*readsimpconf1)
			var m ModelConfParams
			if amino {
				patternval, _ = gophy.PreparePatternVecsProt(t, patternsint, seqs)
				m = ModelConfParams{t: t, amino: amino, aa: *aa}
			} else {
				patternval, _ = gophy.PreparePatternVecs(t, patternsint, seqs)
				m = ModelConfParams{t: t, amino: amino}
			}
			curmodels, curnodemodels, modelnodes, curparams, _ = readModelConf(m)
			curmodels[0].SetBaseFreqs(M.BF) // inherit base model params
			curmodels[0].R = M.R            // for aminos, this is already the same
			curmodels[0].SetupQGTR()
			lm := gophy.PCalcLogLikePatternsMul(t, curmodels, curnodemodels, patternval, *wks)
			fmt.Println("starting lnL:", lm)
			if amino {
				gophy.OptimizeAACompSharedRMNLOPT(t, curmodels, curnodemodels, patternval, false, *wks)
			} else {
				gophy.OptimizeGTRDNACompSharedRM(t, curmodels, curnodemodels, true, patternval, useLog, *wks)
				gophy.OptimizeGTRDNACompSharedRMNL(t, curmodels, curnodemodels, true, patternval, useLog, *wks)
			}
		} else {
			fmt.Println("Calculating likelihood of specified configuration and params")
			t = gophy.ReadTreeFromFile(*readsimpconf1)
			if amino {
				patternval, _ = gophy.PreparePatternVecsProt(t, patternsint, seqs)
			} else {
				patternval, _ = gophy.PreparePatternVecs(t, patternsint, seqs)
			}
			m := ModelConfParams{t: t, modelconf: *readsimpconf2, amino: amino}
			curmodels, curnodemodels, modelnodes, curparams, _ = readModelConf(m)
			//fmt.Println(curmodels, curnodemodels, modelnodes, curparams)
			lm := gophy.PCalcLogLikePatternsMul(t, curmodels, curnodemodels, patternval, *wks)
			currentaic = icfun(lm, curparams, nsites)
			fmt.Println("lm: ", lm, "bic:", currentaic)
			if M.NumStates == 4 {
				gophy.OptimizeGTRDNACompSharedRM(t, curmodels, curnodemodels, true, patternval, useLog, *wks)
				gophy.OptimizeGTRDNACompSharedRMNL(t, curmodels, curnodemodels, true, patternval, useLog, *wks)
			}
			lm = gophy.PCalcLogLikePatternsMul(t, curmodels, curnodemodels, patternval, *wks)
			currentaic = icfun(lm, curparams, nsites)
			fmt.Println("lm: ", lm, "bic:", currentaic)
		}
	} else {
		curmodels, curnodemodels, modelnodes, curparams, currentaic = calcConfiguration(t, &M, amino, *aa, minset,
			patternval, modelparams, l, nsites, useLog, icfun, *wks, *rm, *bl, *gm)
	}

	//at this point we should have the good configuration
	//uncertainty
	// existence, check with and without the
	if *uc1 {
		uncertaintyExist(modelnodes, curmodels, curnodemodels, useLog, t,
			patternval, *wks, curparams, currentaic, nsites, *aicc, amino)
	}
	//uncertainty
	if *uc2 {
		uncertaintyLoc(modelnodes, curmodels, curnodemodels, useLog, t,
			patternval, *wks, curparams, currentaic, nsites, *aicc, amino)
	}
	moveMarksToLabels(t.Rt)
	fmt.Fprintln(os.Stderr, currentaic, curparams)
	fmt.Fprintln(os.Stderr, "Final models")
	fmt.Fprintln(os.Stderr, "-------")
	for i, j := range curmodels {
		fmt.Fprintln(os.Stderr, i, j.BF)
		if !amino && *rm {
			printRM(j.R)
		}
	}
	//make a nexus so it is easier to read in the fig-tree
	fmt.Println(t.Rt.Newick(true) + ";")
	writeFigTreeNexus(curmodels, curnodemodels, t, *tfn+".gophy.results.tre")
}

type ModelConfParams struct {
	t         *gophy.Tree
	modelconf string
	amino     bool
	aa        string
}

func readModelConf(m ModelConfParams) ([]*gophy.DiscreteModel, map[*gophy.Node]int,
	[]*gophy.Node, float64, float64) {
	curparams := 0.
	currentaic := 0.
	allmodels := make(map[int][]*gophy.Node)
	curmodels := make([]*gophy.DiscreteModel, 0)
	curnodemodels := make(map[*gophy.Node]int)
	modelnodes := make([]*gophy.Node, 0)
	if m.modelconf != "" {
		file, err := os.Open(m.modelconf)
		if err != nil {
			log.Fatal(err)
		}
		defer file.Close()
		modelparams := make([]float64, 0)
		reader := bufio.NewReader(file)
		for {
			st, err := reader.ReadString('\n')
			if strings.Contains(st, "bf ") {
				strs := strings.Split(strings.TrimSpace(st), " ")
				//bfn, _ := strconv.Atoi(strs[1][:len(strs[1])-1])
				bfp := make([]float64, 0)
				for _, i := range strs[2:] {
					fp, _ := strconv.ParseFloat(i, 64)
					bfp = append(bfp, fp)
				}
				fmt.Println(bfp)
				var M gophy.DiscreteModel
				if m.amino {
					newmodel := gophy.NewProteinModel()
					newmodel.SetRateMatrixJTT()
					newmodel.M.SetBaseFreqs(bfp)
					newmodel.M.SetupQGTR()
					M = newmodel.M
				} else {
					newmodel := gophy.NewDNAModel()
					newmodel.M.SetBaseFreqs(bfp)
					newmodel.M.SetRateMatrix(modelparams)
					newmodel.M.SetupQGTR()
					M = newmodel.M
				}
				curmodels = append(curmodels, &M)
			}
			if strings.Contains(st, "final: ") {
				strs := strings.Split(strings.TrimSpace(st), " ")
				curparams, _ = strconv.ParseFloat(strs[2], 64)
				currentaic, _ = strconv.ParseFloat(strs[1], 64)
			}
			if strings.Contains(st, "rm: ") {
				strs := strings.Split(strings.TrimSpace(st), " ")
				for _, i := range strs[1:] {
					fp, _ := strconv.ParseFloat(i, 64)
					modelparams = append(modelparams, fp)
				}
			}
			if err == io.EOF {
				break
			}
			if err != nil {
				log.Printf("read %v bytes: %v", st, err)
				break
			}
		}
		for _, i := range m.t.Pre {
			if len(i.Chs) == 0 {
				continue
			}
			nn, err := strconv.Atoi(i.Nam)
			if err != nil {
				fmt.Println("something wrong with the internal node name:", i)
				os.Exit(0)
			}
			curnodemodels[i] = nn
			i.IData["shift"] = nn
			for _, j := range i.Chs {
				if len(j.Chs) == 0 {
					curnodemodels[j] = nn
					j.IData["shift"] = nn
				}
			}
			if i.Par != nil {
				if i.Nam != i.Par.Nam {
					modelnodes = append(modelnodes, i)
				}
			}
		}
	} else {
		for _, i := range m.t.Pre {
			if len(i.Chs) == 0 { // tip
				continue
			}
			nn, err := strconv.Atoi(i.Nam) // node name (here it is model label)
			if err != nil {
				fmt.Println("something wrong with the internal node name:", i)
				os.Exit(0)
			}
			_, ok := allmodels[nn]
			if !ok { // model not in modelset yet
				var M gophy.DiscreteModel
				if m.amino {
					newmodel := gophy.NewProteinModel()
					if m.aa == "JTT" {
						newmodel.SetRateMatrixJTT()
					} else if m.aa == "WAG" {
						newmodel.SetRateMatrixWAG()
					} else if m.aa == "LG" {
						newmodel.SetRateMatrixLG()
					}
					newmodel.M.SetModelBF() // use models freqs for ease, can change later
					newmodel.M.SetupQGTR()
					M = newmodel.M
				} else {
					newmodel := gophy.NewDNAModel()
					bfs := make([]float64, 4)
					for i := range bfs {
						bfs[i] = 1. / 4
					}
					newmodel.M.SetBaseFreqs(bfs) // equal
					modelparams := make([]float64, 5)
					for i := range modelparams {
						modelparams[i] = 1.0
					}
					newmodel.M.SetRateMatrix(modelparams) // equal rates
					newmodel.M.SetupQGTR()
					M = newmodel.M
				}
				allmodels[nn] = append(allmodels[nn], i)
				curmodels = append(curmodels, &M)
			}

			curnodemodels[i] = nn
			i.IData["shift"] = nn
			for _, j := range i.Chs {
				if len(j.Chs) == 0 {
					curnodemodels[j] = nn
					j.IData["shift"] = nn
				}
			}
			if i.Par != nil { // not root
				if i.Nam != i.Par.Nam { // models don't match
					modelnodes = append(modelnodes, i)
				}
			}
		}
	}
	return curmodels, curnodemodels, modelnodes, curparams, currentaic
}

func printRM(R *mat.Dense) {
	for i := 0; i < 4; i++ {
		for j := 0; j < 4; j++ {
			if i < j {
				fmt.Fprint(os.Stderr, R.At(i, j), " ")
			} else {
				fmt.Fprint(os.Stderr, "- ")
			}
		}
		fmt.Fprint(os.Stderr, "\n")
	}
}

func initOptimization(t *gophy.Tree, M *gophy.DiscreteModel, patternval []float64,
	wks int, aa bool, modelparams *[]float64, bl bool, gamma bool) (float64, bool) {
	l := gophy.PCalcLikePatterns(t, M, patternval, wks)
	fmt.Fprintln(os.Stderr, "lnL:", l)
	useLog := false //slower
	if math.IsInf(l, -1) {
		fmt.Println("initial value problem, switch to big")
		useLog = true
		l = gophy.PCalcLogLikePatterns(t, M, patternval, wks)
		fmt.Fprintln(os.Stderr, "lnL:", l)
	}
	if !aa {
		*modelparams = gophy.OptimizeGTRDNA(t, M, patternval, useLog, wks)
	}
	gophy.OptimizeBF(t, M, patternval, useLog, wks)
	if bl {
		if gamma {
			gophy.OptimizeGamma(t, M, patternval, useLog, wks)
			gophy.OptimizeBF(t, M, patternval, useLog, wks)
			gophy.OptimizeGTRDNA(t, M, patternval, useLog, wks)
			gophy.OptimizeGammaBLSNL(t, M, patternval, wks)
			//gophy.OptimizeGammaAndBL(t, M, patternval, useLog, wks)
		} else {
			if useLog {
				gophy.OptimizeBLS(t, M, patternval, wks)
			} else {
				gophy.OptimizeBLNR(t, M, patternval, wks)
			}
		}
	} else {
		if gamma {
			gophy.OptimizeGamma(t, M, patternval, useLog, wks)
			gophy.OptimizeGTRDNA(t, M, patternval, useLog, wks)
			gophy.OptimizeGamma(t, M, patternval, useLog, wks)
			gophy.OptimizeBF(t, M, patternval, useLog, wks)
			gophy.OptimizeGTRDNA(t, M, patternval, useLog, wks)
			gophy.OptimizeBF(t, M, patternval, useLog, wks)
			gophy.OptimizeGTRDNA(t, M, patternval, useLog, wks)
			gophy.OptimizeBF(t, M, patternval, useLog, wks)
			gophy.OptimizeGamma(t, M, patternval, useLog, wks)
			gophy.OptimizeGTRDNA(t, M, patternval, useLog, wks)
			gophy.OptimizeBF(t, M, patternval, useLog, wks)
			gophy.OptimizeGamma(t, M, patternval, useLog, wks)
		}

	}

	if useLog {
		if gamma {
			l = gophy.PCalcLogLikePatternsGamma(t, M, patternval, wks)
		} else {
			l = gophy.PCalcLogLikePatterns(t, M, patternval, wks)
		}
	} else {
		if gamma {
			l = gophy.PCalcLikePatternsGamma(t, M, patternval, wks)
		} else {
			l = gophy.PCalcLikePatterns(t, M, patternval, wks)
		}
	}

	fmt.Fprintln(os.Stderr, "lnL:", l)
	return l, useLog
}

func calcConfiguration(t *gophy.Tree, M *gophy.DiscreteModel, aa bool, aam string, minset int,
	patternval []float64, modelparams []float64, l float64,
	nsites int, useLog bool, icfun func(float64, float64, int) (x float64),
	wks int, rm bool, bl bool, gamma bool) ([]*gophy.DiscreteModel, map[*gophy.Node]int,
	[]*gophy.Node, float64, float64) {

	//for each node in the tree if the number of tips is > minset
	allmodels := make([]*gophy.DiscreteModel, 1) //number of clades that have enough taxa and aren't the root
	modelmap := make(map[*gophy.Node]int)
	allmodels[0] = M
	count := 1
	for _, i := range t.Post {
		if i == t.Rt {
			continue
		}
		if len(i.GetTips()) > minset {
			if aa {
				y := gophy.NewProteinModel()
				y.M.SetBaseFreqs(M.BF)
				if aam == "JTT" {
					y.SetRateMatrixJTT()
				} else if aam == "LG" {
					y.SetRateMatrixLG()
				} else if aam == "WAG" {
					y.SetRateMatrixWAG()
				}
				y.M.SetupQGTR()
				if gamma {
					y.M.GammaNCats = M.GammaNCats
					y.M.GammaCats = M.GammaCats
					y.M.GammaAlpha = M.GammaAlpha
				}
				modelmap[i] = count
				allmodels = append(allmodels, &y.M)
			} else {
				y := gophy.NewDNAModel()
				y.M.SetBaseFreqs(M.BF)
				y.M.SetRateMatrix(modelparams)
				y.M.SetupQGTR()
				if gamma {
					y.M.GammaNCats = M.GammaNCats
					y.M.GammaCats = M.GammaCats
					y.M.GammaAlpha = M.GammaAlpha
				}
				modelmap[i] = count
				allmodels = append(allmodels, &y.M)
			}
			count++
		}
	}
	fmt.Fprintln(os.Stderr, "num of models:", count-1)

	nodemodels := make(map[*gophy.Node]int)
	numbaseparams := ((2 * float64(len(t.Tips))) - 3.) + float64(M.NumStates) - 1
	if !aa {
		numbaseparams += 5. //and GTR
	}
	if gamma {
		numbaseparams += 1.
	}
	nodevalues := make(map[*gophy.Node]float64) // key node, value aicc
	saic := icfun(l, numbaseparams, nsites)
	currentaic := saic
	//start with the sorted set
	fmt.Fprintln(os.Stderr, "starting IC:", saic)
	//fmt.Fprintln(os.Stderr, useLog)
	cur := 1
	for _, i := range t.Post {
		if i == t.Rt {
			continue
		}
		if getVisible(i) > minset {
			fmt.Fprint(os.Stderr, "\rOn "+strconv.Itoa(cur)+"/"+strconv.Itoa(count-1), " ")
			//fmt.Fprint(os.Stderr, i.Newick(false))
			y := allmodels[modelmap[i]]
			models := []*gophy.DiscreteModel{allmodels[0], y}
			if aa || !rm {
				gophy.OptimizeBFSubCladeNLOPT(t, i, false, y, patternval, useLog, wks)
			} else {
				gophy.OptimizeBFDNARMSubClade(t, i, false, y, patternval, wks)
			}

			for _, j := range t.Post {
				nodemodels[j] = 0
			}
			for _, j := range i.PreorderArray() {
				nodemodels[j] = 1
			}
			lm := 1.0
			if useLog {
				lm = gophy.PCalcLogLikePatternsMul(t, models, nodemodels, patternval, wks)
			} else {
				lm = gophy.PCalcLikePatternsMul(t, models, nodemodels, patternval, wks)
			}
			newparamnums := numbaseparams + float64(allmodels[0].NumStates-1.)
			if rm {
				newparamnums += (((float64(allmodels[0].NumStates) * float64(allmodels[0].NumStates)) - float64(allmodels[0].NumStates)) / 2.) - 1.
			}
			//fmt.Println(lm)
			naic := icfun(lm, newparamnums, nsites)
			//fmt.Println(i.Newick(false))
			//fmt.Println("bic: ", lm, naic, y.BF)
			//os.Exit(0)
			nodevalues[i] = naic
			cur++
		}
	}
	keys := sortAicMap(nodevalues)
	curmodels := []*gophy.DiscreteModel{allmodels[0]}
	curnodemodels := make(map[*gophy.Node]int)
	//get the final model configuration
	for _, j := range t.Post {
		curnodemodels[j] = 0
	}
	fmt.Fprintln(os.Stderr, "\rfinal est")
	curparams := numbaseparams
	modelnodes := []*gophy.Node{}
	cur = 1
	for _, k := range keys {
		fmt.Fprintln(os.Stderr, "On "+strconv.Itoa(cur)+"/"+strconv.Itoa(count-1)+" ("+fmt.Sprintf("%f", k.Value)+")")
		//skip if going to take over root model
		if getVisible(t.Rt)-getVisible(k.Key) <= minset {
			continue
		}
		if getVisible(k.Key) > minset {
			//fmt.Println(k.Key, k.Value)
			//shortcut, if k.Value is worse than the best by > 10 don't consider it
			if k.Value-saic > 35 { // 10 is arbitrary. pick another measure
				fmt.Fprintln(os.Stderr, "   ", k.Value)
				cur++
				continue
			}
			//fmt.Println(k.Value, k.Key)
			y := allmodels[modelmap[k.Key]]
			testmodels, testnodemodels, modelint := getNodeModels(curmodels, curnodemodels, y, k.Key)
			//need to be able to send better starting points so it doesn't take as long
			lm := 1.0
			if aa {
				gophy.OptimizeAACompSharedRMNLOPT(t, testmodels, testnodemodels, patternval, false, wks)
			} else {
				if rm {
					gophy.OptimizeGTRBPDNAMul(t, testmodels, testnodemodels, true, patternval, useLog, wks)
				} else {
					gophy.OptimizeGTRDNACompSharedRM(t, testmodels, testnodemodels, true, patternval, useLog, wks)
					gophy.OptimizeGTRDNACompSharedRMNL(t, testmodels, testnodemodels, true, patternval, useLog, wks)
				}
				if bl {
					if gamma {
						gophy.OptimizeSharedGammaAndBLMult(t, testmodels, testnodemodels, patternval, useLog, wks)
					} else {
						gophy.OptimizeBLNRMult(t, testmodels, testnodemodels, patternval, wks)
					}
				}
			}
			if useLog {
				lm = gophy.PCalcLogLikePatternsMul(t, testmodels, testnodemodels, patternval, wks)
			} else {
				lm = gophy.PCalcLikePatternsMul(t, testmodels, testnodemodels, patternval, wks)
			}
			addparams := float64(testmodels[0].NumStates) - 1.
			if rm {
				addparams += (((float64(allmodels[0].NumStates) * float64(allmodels[0].NumStates)) - float64(allmodels[0].NumStates)) / 2.) - 1.
			}
			naic := icfun(lm, curparams+addparams, nsites)
			fmt.Fprintln(os.Stderr, " -- ", naic, lm, currentaic)
			if naic < currentaic {
				fmt.Fprintln(os.Stderr, naic, "<", currentaic, curparams+addparams)
				curparams += addparams
				curmodels = testmodels
				curnodemodels = testnodemodels
				modelnodes = append(modelnodes, k.Key)
				if bl {
					if gamma {
						gophy.OptimizeSharedGammaAndBLMult(t, curmodels, curnodemodels, patternval, true, wks)
					} else {
						gophy.OptimizeBLNRMult(t, curmodels, curnodemodels, patternval, wks)
					}
					lm = gophy.PCalcLogLikePatternsMul(t, curmodels, curnodemodels, patternval, wks)
					naic = icfun(lm, curparams, nsites)
					fmt.Fprintln(os.Stderr, " (final bl)", lm, naic)
				}
				currentaic = naic
			} else {
				unmarkModel(k.Key, modelint)
			}
		}
		cur++
	}
	return curmodels, curnodemodels, modelnodes, curparams, currentaic
}

func uncertaintyLoc(modelnodes []*gophy.Node, curmodels []*gophy.DiscreteModel,
	curnodemodels map[*gophy.Node]int, useLog bool, t *gophy.Tree, patternval []float64,
	wks int, curparams float64, currentaic float64, nsites int, useaicc bool, aa bool) { // location
	var icfun func(float64, float64, int) (x float64)
	icfun = gophy.CalcBIC
	if useaicc {
		icfun = gophy.CalcAICC
	}
	mainaic := 1.0 //same as math.Exp((currentaic - currentaic) / 2)
	for _, i := range modelnodes {
		fmt.Fprintln(os.Stderr, i, i.IData["shift"])
		testnodes := []*gophy.Node{}
		if _, ok := i.Par.IData["shift"]; !ok {
			testnodes = append(testnodes, i.Par)
		}
		for _, j := range i.Chs {
			if j.IData["shift"] == i.IData["shift"] {
				testnodes = append(testnodes, j)
			}
		}
		tevals := []float64{mainaic}
		for _, ts := range testnodes {
			//unmark just the modelnodes
			testnodemodels := make(map[*gophy.Node]int)
			testmodels := []*gophy.DiscreteModel{}
			for _, j := range curmodels {
				testmodels = append(testmodels, j.DeepCopyDiscreteModel())
			}
			for k, j := range curnodemodels {
				testnodemodels[k] = j
			}
			if ts != i.Par {
				if _, ok := i.Par.IData["shift"]; !ok {
					testnodemodels[i] = 0
				} else {
					testnodemodels[i] = testnodemodels[i.Par]
				}
				//get other child
				for _, k := range ts.GetSib().PreorderArray() {
					testnodemodels[k] = 0
				}
			} else { // the node is the parent
				testnodemodels[ts] = i.IData["shift"]
				for _, k := range ts.PreorderArray() {
					if _, ok := k.IData["shift"]; !ok {
						testnodemodels[k] = i.IData["shift"]
					}
				}
			}
			//get aic, could reestimate the model but not right now, if this happens it needs to be copied so that we can
			//   go back to the original as well
			tlm := 0.0
			if aa {
				gophy.OptimizeAACompSharedRMSingleModel(t, testmodels,
					testnodemodels, true, i.IData["shift"], patternval, useLog, wks)
			} else {
				gophy.OptimizeGTRDNACompSharedRMSingleModel(t, testmodels,
					testnodemodels, true, i.IData["shift"], patternval, useLog, wks)
			}
			if useLog {
				tlm = gophy.PCalcLogLikePatternsMul(t, testmodels, testnodemodels, patternval, wks)
			} else {
				tlm = gophy.PCalcLikePatternsMul(t, testmodels, testnodemodels, patternval, wks)
			}
			taic := icfun(tlm, curparams, nsites)
			tstat := math.Exp((currentaic - taic) / 2)
			tevals = append(tevals, tstat)
		}
		waic := mainaic / floats.Sum(tevals)
		i.FData["uncloc"] = waic
		fmt.Fprintln(os.Stderr, currentaic, waic)
	}
}

func uncertaintyExist(modelnodes []*gophy.Node, curmodels []*gophy.DiscreteModel,
	curnodemodels map[*gophy.Node]int, useLog bool, t *gophy.Tree, patternval []float64,
	wks int, curparams float64, currentaic float64, nsites int, useaicc bool, aa bool) {
	var icfun func(float64, float64, int) (x float64)
	icfun = gophy.CalcBIC
	if useaicc {
		icfun = gophy.CalcAICC
	}
	mainaic := 1.0 //same as math.Exp((currentaic - currentaic) / 2)
	for _, i := range modelnodes {
		fmt.Fprintln(os.Stderr, i, i.IData["shift"])
		//unmark just the modelnodes
		testnodemodels := make(map[*gophy.Node]int)
		testmodels := []*gophy.DiscreteModel{}
		for _, j := range curmodels {
			testmodels = append(testmodels, j.DeepCopyDiscreteModel())
		}
		for k, j := range curnodemodels {
			testnodemodels[k] = j
		}
		for _, k := range i.PreorderArray() {
			if i.IData["shift"] == k.IData["shift"] {
				testnodemodels[k] = 0
			}
		}
		//get aic, could reestimate the model but not right now, if this happens it needs to be copied so that we can
		//   just reestimate one model
		//   go back to the original as well
		tlm := 0.0
		if aa {
			gophy.OptimizeAACompSharedRMSingleModel(t, testmodels,
				testnodemodels, true, 0, patternval, useLog, wks)
		} else {
			//gophy.OptimizeGTRDNACompSharedRMSingleModel(t, testmodels,
			//	testnodemodels, true, 0, patternval, useLog, wks)
			gophy.OptimizeGTRDNACompSharedRM(t, testmodels, testnodemodels, true, patternval, useLog, wks)
			gophy.OptimizeGTRDNACompSharedRMNL(t, testmodels, testnodemodels, true, patternval, useLog, wks)

		}
		if useLog {
			tlm = gophy.PCalcLogLikePatternsMul(t, testmodels, testnodemodels, patternval, wks)
		} else {
			tlm = gophy.PCalcLikePatternsMul(t, testmodels, testnodemodels, patternval, wks)
		}
		addparams := float64(testmodels[0].NumStates) - 1.
		taic := icfun(tlm, curparams-addparams, nsites)
		tstat := math.Exp((currentaic - taic) / 2)
		i.FData["uncex"] = mainaic / (mainaic + tstat)
		fmt.Fprintln(os.Stderr, currentaic, taic, "(", tlm, ")", i.FData["uncex"])
		fmt.Fprintln(os.Stderr, "exists:", i.IData["shift"], i.FData["uncex"])
	}
}

func writeFigTreeNexus(curmodels []*gophy.DiscreteModel, curnodemodels map[*gophy.Node]int, t *gophy.Tree,
	outfile string) {
	f, err := os.Create(outfile)
	if err != nil {
		panic(err)
	}
	w := bufio.NewWriter(f)
	w.WriteString("#NEXUS\nbegin trees;\ntree t = ")
	delim := ","
	for _, i := range t.Post {
		i.IData["model"] = curnodemodels[i]
		i.SData["modelpar"] = "{" + strings.Trim(strings.Join(strings.Fields(fmt.Sprint(curmodels[curnodemodels[i]].BF)), delim), "[]") + "}"
	}
	X := nexusRecur(t.Rt, true)
	w.WriteString(X + ";\nend;\n")
	w.Flush()
}

func nexusRecur(n *gophy.Node, bl bool) (ret string) {
	var buffer bytes.Buffer
	for in, cn := range n.Chs {
		if in == 0 {
			buffer.WriteString("(")
		}
		buffer.WriteString(nexusRecur(cn, bl))
		if bl {
			s := strconv.FormatFloat(cn.Len, 'f', -1, 64)
			buffer.WriteString(":")
			buffer.WriteString(s)
		}
		if in == len(n.Chs)-1 {
			buffer.WriteString(")")
		} else {
			buffer.WriteString(",")
		}
	}
	//CHANGE THIS PART
	if len(n.Chs) == 0 {
		buffer.WriteString(n.Nam)
	}
	buffer.WriteString("[&" + "model" + "=" + strconv.Itoa(n.IData["model"]))
	buffer.WriteString(",modelpar=" + n.SData["modelpar"])
	if _, ok := n.FData["uncloc"]; ok {
		buffer.WriteString(",uncloc=" + strconv.FormatFloat(n.FData["uncloc"], 'f', -1, 64))
	}
	if _, ok := n.FData["uncex"]; ok {
		buffer.WriteString(",uncex=" + strconv.FormatFloat(n.FData["uncex"], 'f', -1, 64))
	}
	buffer.WriteString("]")
	//END CHANGE
	ret = buffer.String()
	return
}
