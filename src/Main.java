import java.util.ArrayList;
import java.util.Collections;

/**
 * <n> This class is to test if the sample median is an unbiased estimator of the population mean. </n>
 *
 * It finds the sample median is an unbiased estimator of the population mean.
 */
public class Main {
    public static void main(String[] args) {
        // How large the population should be
        int POPULATION_SIZE = 1000;

        /* How large the sample is, as a percentage of the population size. */
        double PCT_SAMPLE_SIZE = 0.05;

        // Computes the sample size
        int SAMPLE_SIZE = (int) ((double)(POPULATION_SIZE) * PCT_SAMPLE_SIZE);

        // How many trials to run
        double NUM_TRIALS = 100000;

        // Counts the number of times the sample median is larger than the population mean
        double COUNT_SAMPLE_MEDIAN_GREATER = 0;

        // Compare if the sample median is higher than the population mean
        // If yes: increment the first count variable
        for(int i = 0; i < NUM_TRIALS; i++){
            if(runTrial(POPULATION_SIZE, SAMPLE_SIZE)){
                COUNT_SAMPLE_MEDIAN_GREATER++;
            }
        }

        // Determines the percentage of time the sample median is larger than the population mean
        double PCT_SAMPLE_MEDIAN_GREATER = COUNT_SAMPLE_MEDIAN_GREATER / NUM_TRIALS;

        // Create a confidence interval for P(sample median > population mean).
        ArrayList<Double> confidenceInterval = getConfidenceInterval(PCT_SAMPLE_MEDIAN_GREATER, (int) NUM_TRIALS);

        // Check if this confidence interval contains P = 50%.
        // If yes: the sample median is an unbiased estimator of the population mean
        // Else: the sample median is a biased estimator of the population mean

        String s = "The sample median is ";
        if(!containsFiftyPercent(confidenceInterval)){
            s += " NOT ";
        }
        s += "an unbiased estimator of the population mean.";

        System.out.println(s);
        System.out.println("Confidence interval (n = " + (int) NUM_TRIALS + "): " + confidenceInterval);
    }

    /**
     * Checks if the given confidence interval for a proportion contains the value X = 50%.
     *
     * @param confidenceInterval the confidence interval to test
     * @return whether the value X = 50% is included in the bounds.
     * @throws RuntimeException if the confidence interval array is empty.
     */
    public static boolean containsFiftyPercent(ArrayList<Double> confidenceInterval){
        if(confidenceInterval.isEmpty()){
            throw new RuntimeException("The given array is empty!");
        }

        double lowerBound = confidenceInterval.get(0);
        double upperBound = confidenceInterval.get(1);

        return lowerBound <= 0.5 && 0.5 <= upperBound;
    }

    /**
     * Returns the 95% confidence interval for the given sample proportion.
     * @param sampleProportion the sample proportion
     * @param sampleSize the number of sample observations
     * @return the 95% CI for the population proportion.
     * @throws RuntimeException if either input is zero or less
     */
    private static ArrayList<Double> getConfidenceInterval(double sampleProportion, int sampleSize){
        if(sampleProportion <= 0 || sampleSize <= 0){
            throw new RuntimeException("At least one of the inputs are non-positive, cannot run calculations.");
        }

        ArrayList<Double> output = new ArrayList<>();

        double marginOfError = getMarginOfError(sampleProportion, sampleSize);

        double lowerBound = sampleProportion - marginOfError;
        double upperBound = sampleProportion + marginOfError;

        output.add(lowerBound); output.add(upperBound);

        return output;
    }

    /**
     * Calculates the margin of error for a sample of proportions.
     * @param sampleMean the mean of the sample
     * @param sampleSize the number of observations of the sample
     * @return the margin of error for the true population proportion estimate, at 95% confidence.
     * @throws RuntimeException if the sample size is less than 1
     */
    private static double getMarginOfError(double sampleMean, int sampleSize){
        if(sampleSize <= 0){
            throw new RuntimeException("Sample size cannot be less than 1!");
        }

        return 1.96 * getStandardError(sampleMean, sampleSize);
    }

    /**
     * Calculates the standard error for a proportion.
     * @param sampleProportion the mean of the sample
     * @param sampleSize the number of sample observations
     * @return the standard error of the sample
     * @throws RuntimeException if the sample size is less than 1
     */
    private static double getStandardError(double sampleProportion, int sampleSize){
        if(sampleSize <= 0){
            throw new RuntimeException("Sample size cannot be less than 1!");
        }

        return Math.sqrt((sampleProportion * (1 - sampleProportion))/(double) sampleSize);
    }

    /**
     * Conducts a single trial.
     *
     * <n> A trial consists of creating a random population,
     * seleecting a sample, and
     * computing both a sample median and a population mean.
     * The sample median is compared to the population mean. </n>
     *
     * @param POPULATION_SIZE the number of elements in the population
     * @param SAMPLE_SIZE the number of elements in the sample
     *
     * @return whether the sample median is higher than the population mean
     * @throws RuntimeException if either input is less than 1, or if the population size is smaller than the sample.
     */
    public static boolean runTrial(int POPULATION_SIZE, int SAMPLE_SIZE){
        if(SAMPLE_SIZE <= 0 || POPULATION_SIZE < SAMPLE_SIZE){
            throw new RuntimeException("At least one of the inputs is less than 1, or the population size is smaller than the sample.");
        }

        // Generate a population
        ArrayList<Double> population = initPopulation(POPULATION_SIZE);

        // Grab a sample
        ArrayList<Double> sample = getSample(population, SAMPLE_SIZE);

        // Calculate the sample median
        double sampleMedian = getMedian(sample);

        // Calculate the population mean
        double popMean = getAverage(population);

        return sampleMedian > popMean;
    }

    /**
     * Returns the median of the given array.
     * @throws RuntimeException if the given array is empty.
     */
    private static double getMedian(ArrayList<Double> arr) {
        if(arr.isEmpty()){
            throw new RuntimeException("The given array is empty!");
        }
        
        Collections.sort(arr);
        // sanitizes the input before operations

        double size = arr.size();
        double half_size = size / 2;

        int HALF_SIZE_LOWER_BOUND = (int)(half_size) - 1;
        int HALF_SIZE_UPPER_BOUND = (int)(half_size);

        if(size % 2 == 0)
        {
            double a1 = arr.get(HALF_SIZE_LOWER_BOUND);
            double a2 = arr.get(HALF_SIZE_UPPER_BOUND);

            return getAverage(a1, a2);
        } // IF even? find the middle two values + average them!
        else
        {
            return arr.get((int)(half_size));
        } // ELSE? returns the value in the middle of the ArrayList
    }

    public static double getAverage(double d1, double d2){
        return 0.5 * (d1 + d2);
    }

    /**
     * Returns a random double, min- and max-inclusive.
     * @param min the lower bound of the output
     * @param max the upper bound of the output
     * @return a random double in the range
     */
    private static double getRandDouble(double min, double max)
    {
        return Math.random() * (max - min + 1) + min;
    }

    /**
     * Returns an array with uniformly distributed random doubles across X = [0.0, 100.0].
     * @param populationSize the number of elements to store in the array.
     * @throws RuntimeException if the population size is less than 1
     */
    private static ArrayList<Double> initPopulation(int populationSize) {
        if(populationSize <= 0){
            throw new RuntimeException("Population size cannot be less than 1!");
        }

        ArrayList<Double> output = new ArrayList<>();

        for(int i = 0; i < populationSize; i++){
            double randRating = getRandDouble(0, 100);
            output.add(randRating);
        }

        return output;
    }

    /**
     * Returns the mean of the given array.
     * @throws RuntimeException if the array is empty.
     */
    private static double getAverage(ArrayList<Double> arr) {
        if(arr.isEmpty()){
            throw new RuntimeException("The given array is empty!");
        }
        return getSum(arr) / arr.size();
    }

    /**
     * Returns the sum of the given array.
     * @throws RuntimeException if the given array is empty.
     */
    private static double getSum(ArrayList<Double> arr){
        if(arr.isEmpty()){
            throw new RuntimeException("The given array is empty!");
        }
        
        double sum = 0;
        for(Double d : arr){
            sum += d;
        }
        return sum;
    }

    /**
     * Returns a sample out of the given population.
     * @param population the population to sample from
     * @param sampleSize the number of elements the sample should have.
     * @throws RuntimeException if the sample size is larger than the population array
     */
    private static ArrayList<Double> getSample(ArrayList<Double> population, int sampleSize) {
        if(sampleSize > population.size()){
            throw new RuntimeException("The sample size can't be larger than the population!");
        }

        ArrayList<Double> output = new ArrayList<>();
        Collections.shuffle(population);

        for(int i = 0; i < sampleSize; i++){
            output.add(population.get(i));
        }

        return output;
    }
}