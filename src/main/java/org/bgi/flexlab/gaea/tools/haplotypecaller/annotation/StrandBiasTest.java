package org.bgi.flexlab.gaea.tools.haplotypecaller.annotation;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.StringTokenizer;

import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.tools.haplotypecaller.ReadLikelihoods;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;
import org.bgi.flexlab.gaea.util.Utils;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;

public abstract class StrandBiasTest extends InfoFieldAnnotation {

    protected static final int ARRAY_DIM = 2;
    protected static final int ARRAY_SIZE = ARRAY_DIM * ARRAY_DIM;

    @Override
    //template method for calculating strand bias annotations using the three different methods
    public Map<String, Object> annotate(final ChromosomeInformationShare ref,
                                        final VariantContext vc,
                                        final ReadLikelihoods<Allele> likelihoods) {
        Utils.nonNull(vc);
        if ( !vc.isVariant() ) {
            return Collections.emptyMap();
        }

        // if the genotype and strand bias are provided, calculate the annotation from the Genotype (GT) field
        if ( vc.hasGenotypes() ) {
            for (final Genotype g : vc.getGenotypes()) {
                if (g.hasAnyAttribute(GaeaVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY)) {
                    return calculateAnnotationFromGTfield(vc.getGenotypes());
                }
            }
        }

        if (likelihoods != null) {
            return calculateAnnotationFromLikelihoods(likelihoods, vc);
        }
        return Collections.emptyMap();
    }

    protected abstract Map<String, Object> calculateAnnotationFromGTfield(final GenotypesContext genotypes);

    protected abstract Map<String, Object> calculateAnnotationFromLikelihoods(final ReadLikelihoods<Allele> likelihoods,
                                                                              final VariantContext vc);

    /**
     * Create the contingency table by retrieving the per-sample strand bias annotation and adding them together
     * @param genotypes the genotypes from which to pull out the per-sample strand bias annotation
     * @param minCount minimum threshold for the sample strand bias counts for each ref and alt.
     *                 If both ref and alt counts are above minCount the whole sample strand bias is added to the resulting table
     * @return the table used for several strand bias tests, will be null if none of the genotypes contain the per-sample SB annotation
     */
    protected int[][] getTableFromSamples( final GenotypesContext genotypes, final int minCount ) {
        if( genotypes == null ) {
            return null;
        }

        final int[] sbArray = {0,0,0,0}; // reference-forward-reverse -by- alternate-forward-reverse
        boolean foundData = false;

        for( final Genotype g : genotypes ) {
            if( g.isNoCall() || ! g.hasAnyAttribute(GaeaVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY) ) {
                continue;
            }

            foundData = true;
            int[] data;
            if ( g.getAnyAttribute(GaeaVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY).getClass().equals(String.class)) {
                final String sbbsString = (String)g.getAnyAttribute(GaeaVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY);
                data = encodeSBBS(sbbsString);
            } else if (g.getAnyAttribute(GaeaVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY).getClass().equals(ArrayList.class)) {
                @SuppressWarnings("unchecked")
                final List<Integer> sbbsList = (ArrayList<Integer>) g.getAnyAttribute(GaeaVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY);
                data = encodeSBBS(sbbsList);
            } else {
                throw new UserException("Unexpected " + GaeaVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY + " type");
            }

            if ( passesMinimumThreshold(data, minCount) ) {
                for( int index = 0; index < sbArray.length; index++ ) {
                    sbArray[index] += data[index];
                }
            }
        }

        return foundData ? decodeSBBS(sbArray) : null ;
    }

    /**
     Allocate and fill a 2x2 strand contingency table.  In the end, it'll look something like this:
     *             fw      rc
     *   allele1   #       #
     *   allele2   #       #
     * @return a 2x2 contingency table
     */
    public static int[][] getContingencyTable( final ReadLikelihoods<Allele> likelihoods,
                                               final VariantContext vc,
                                               final int minCount) {
        return getContingencyTable(likelihoods, vc, minCount, likelihoods.samples());
    }

    /**
     Allocate and fill a 2x2 strand contingency table.  In the end, it'll look something like this:
     *             fw      rc
     *   allele1   #       #
     *   allele2   #       #
     * @return a 2x2 contingency table
     */
    public static int[][] getContingencyTable( final ReadLikelihoods<Allele> likelihoods,
                                               final VariantContext vc,
                                               final int minCount,
                                               final Collection<String> samples) {
        if( likelihoods == null || vc == null) {
            return null;
        }

        final Allele ref = vc.getReference();
        final List<Allele> allAlts = vc.getAlternateAlleles();

        final int[][] table = new int[ARRAY_DIM][ARRAY_DIM];
        for (final String sample : samples) {
            final int[] sampleTable = new int[ARRAY_SIZE];
            likelihoods.bestAlleles(sample).stream()
                    .filter(ba -> ba.isInformative())
                    .forEach(ba -> updateTable(sampleTable, ba.allele, ba.read, ref, allAlts));
            if (passesMinimumThreshold(sampleTable, minCount)) {
                copyToMainTable(sampleTable, table);
            }
        }

        return table;
    }

    /**
     * Helper method to copy the per-sample table to the main table
     *
     * @param perSampleTable   per-sample table (single dimension)
     * @param mainTable        main table (two dimensions)
     */
    private static void copyToMainTable(final int[] perSampleTable, final int[][] mainTable) {
        mainTable[0][0] += perSampleTable[0];
        mainTable[0][1] += perSampleTable[1];
        mainTable[1][0] += perSampleTable[2];
        mainTable[1][1] += perSampleTable[3];
    }

    private static void updateTable(final int[] table, final Allele allele, final GaeaSamRecord read, final Allele ref, final List<Allele> allAlts) {
        final boolean matchesRef = allele.equals(ref, true);
        final boolean matchesAnyAlt = allAlts.contains(allele);

        if ( matchesRef || matchesAnyAlt ) {
            final int offset = matchesRef ? 0 : ARRAY_DIM;

            // a normal read with an actual strand
            final boolean isFW = !read.getReadNegativeStrandFlag();
            table[offset + (isFW ? 0 : 1)]++;
        }
    }

    /**
     * Does this strand data array pass the minimum threshold for inclusion?
     *
     * @param data  the array
     * @param minCount The minimum threshold of counts in the array
     * @return true if it passes the minimum threshold, false otherwise
     */
    protected static boolean passesMinimumThreshold(final int[] data, final int minCount) {
        // the ref and alt totals must be greater than MIN_COUNT
        return data[0] + data[1] + data[2] + data[3] > minCount;
    }

    /**
     * Helper function to parse the genotype annotation into the SB annotation array
     * @param string the string that is returned by genotype.getAnnotation("SB")
     * @return the array used by the per-sample Strand Bias annotation
     */
    private static int[] encodeSBBS( final String string ) {
        final int[] array = new int[ARRAY_SIZE];
        final StringTokenizer tokenizer = new StringTokenizer(string, ",", false);
        for( int index = 0; index < ARRAY_SIZE; index++ ) {
            array[index] = Integer.parseInt(tokenizer.nextToken());
        }
        return array;
    }

    /**
     * Helper function to parse the genotype annotation into the SB annotation array
     * @param arrayList the ArrayList returned from StrandBiasBySample.annotate()
     * @return the array used by the per-sample Strand Bias annotation
     */
    private static int[] encodeSBBS( final List<Integer> arrayList ) {
        final int[] array = new int[ARRAY_SIZE];
        int index = 0;
        for ( Integer item : arrayList ) {
            array[index++] = item;
        }

        return array;
    }

    /**
     * Helper function to turn the  SB annotation array into a contingency table
     * @param array the array used by the per-sample Strand Bias annotation
     * @return the table used by the StrandBiasTest annotation
     */
    private static int[][] decodeSBBS( final int[] array ) {
        if(array.length != ARRAY_SIZE) {
            return null;
        }
        final int[][] table = new int[ARRAY_DIM][ARRAY_DIM];
        table[0][0] = array[0];
        table[0][1] = array[1];
        table[1][0] = array[2];
        table[1][1] = array[3];
        return table;
    }
}
