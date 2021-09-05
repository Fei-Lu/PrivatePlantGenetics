package analysis.wheat.VMap2.build;

public class VMap2BuildGo {

    public VMap2BuildGo() {
//        this.firstBuild();
        this.secondBuild();
    }

    public void secondBuild() {
//        this.depthProfile2();
//        this.taxaValidation();
        this.filterSNPs();
//        this.filterGenes();
    }

    public void filterGenes () {
        new FilterGenes();
    }

    public void filterSNPs () {
        new FilterSNPs();
    }

    public void taxaValidation () {
        new TaxaValidation();
    }

    public void depthProfile2 () {
        new DepthProfile2();
    }

    public void firstBuild () {
        //        this.depthProfile();
//        this.filterByDepth();
//        this.scanLibrary();
//        this.filterMafAndOccurrence();
//        this.indelCheck();
    }

    public void indelCheck () {
        new IndelCheck();
    }

    public void filterMafAndOccurrence () {
        new MafAndOccurrence();
    }

    public void scanLibrary () {
        new ScanLibrary();
    }

    public void filterByDepth () {
        new FilterByDepth();
    }

    public void depthProfile () {
        new DepthProfile();
    }

}
