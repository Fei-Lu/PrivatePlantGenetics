package learning;

public class LearnJar {
    public LearnJar () {
        this.createJar();
        this.viewJar();
        this.extractJar();
        this.updateJar();
        this.modifyManifest();
        this.signAndVerify();
    }

    public void signAndVerify () {
        //learn when needed.
    }

    public void modifyManifest () {
        //run the command line in terminal
        String command = "jar cfm /Users/feilu/Documents/analysisL/learning/jar/learnJar.jar /Users/feilu/Documents/NetBeansProjects/PlantGenetics/PrivatePlantGenetics/src/learning/MANIFEST.MF /Users/feilu/Documents/NetBeansProjects/PlantGenetics/PrivatePlantGenetics/src/learning/LearnJar.java";
    }

    public void updateJar () {
        //run the command line in terminal
        String command = "jar uf /Users/feilu/Documents/analysisL/learning/jar/learnJar.jar /Users/feilu/Documents/NetBeansProjects/PlantGenetics/PrivatePlantGenetics/src/learning/LearnSection.java";

    }

    public void extractJar () {
        //run the command line in terminal
        String command = "jar xf /Users/feilu/Documents/analysisL/learning/jar/learnJar.jar";
    }

    public void viewJar () {
        //run the command line in terminal
        String command = "jar tf /Users/feilu/Documents/analysisL/learning/jar/learnJar.jar";
    }

    public void createJar () {
        //run the command line in terminal
        String command = "jar cvf /Users/feilu/Documents/analysisL/learning/jar/learnJar.jar /Users/feilu/Documents/NetBeansProjects/PlantGenetics/PrivatePlantGenetics/src/learning/LearnJar.java";
    }

    public static void main (String[] args) {
        System.out.println("This is LearJar!");
    }
}
