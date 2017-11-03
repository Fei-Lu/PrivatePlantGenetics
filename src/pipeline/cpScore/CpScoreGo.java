/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pipeline.cpScore;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

/**
 *
 * @author feilu
 */
public class CpScoreGo {
    Options options = new Options();
    HelpFormatter optionFormat = new HelpFormatter();
    String usage = null;
    int kmerLength = 0;
    
    public CpScoreGo (String[] args) {
        usage = this.getUsage();
        this.buildOptions();
        this.outputUsage();
        this.initialize(args);
    }
    
    
    void initialize (String[] args) {
        CommandLineParser parser = new DefaultParser();
        try {
            CommandLine line = parser.parse( options, args);
            String value = null;
            if( line.hasOption( "m" ) ) {
                value = line.getOptionValue("m");
                System.out.println( line.getOptionValue( "l" ) );
            }
        }
        catch(ParseException exp ) {
            System.out.println( "Unexpected exception:" + exp.getMessage() );
        }
    }
    
    void outputUsage () {
        optionFormat.printHelp( usage, options );
    }
    
    void buildOptions () {
        options = new Options();
        options.addOption("m", true, "Analysis mode. Two modes are available, building reference kmer library (b option) and profiling CpScore (p option). e.g. -m b");
        options.addOption("k", true, "Kmer length. Only 32 and 16 are supported. e.g. -k 32");
        options.addOption("i", true, "Input file. e.g -i maizeAGPV4.fa");
        options.addOption("l", true, "Kmer library file. e.g -l maize_32mer.lib");
        options.addOption("o", true, "Output file or directory. e.g -o maize_32mer.lib");
    }
    
    String getUsage () {
        StringBuilder sb = new StringBuilder();
        sb.append("\nThe program CpScoreProfiler.jar is designed to calculate base copy number of a given genome.");
        sb.append("It has 2 seperate analysis modes. The first is to build kmer library from the reference genome. The second is to calculate the CpScore from another non-reference genome.\n\n");
        sb.append("Command line example:\n");
        sb.append("\t1. Build kmer library from reference genome.\n");
        sb.append("\tjava -jar CpScoreProfiler.jar -m b -k 32 -i maizeAGPV4.fa -o maize_32mer.lib\n");
        sb.append("\t2. Calculate CpScore from another non-reference genome\n");
        sb.append("\tjava -jar CpScoreProfiler.jar -m p -k 32 -l maize_32mer.lib -i CML247.fa -o CML247_Cp\n\n");
        return sb.toString();
    }
    
    public static void main (String[] args) { 
        new CpScoreGo (args);
    }
}
