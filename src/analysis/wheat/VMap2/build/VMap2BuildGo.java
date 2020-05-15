package analysis.wheat.VMap2.build;

public class VMap2BuildGo {

    public VMap2BuildGo() {
        //this.depthProfile();
        //this.filterByDepth();
        this.scanLibrary();
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
