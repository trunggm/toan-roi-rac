package im.cleGrasp;

import java.util.Vector;

import com.sun.org.apache.regexp.internal.recompile;

import sun.security.util.DisabledAlgorithmConstraints;
import sun.security.util.Length;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;


public class Compute {
	
	public Compute () {
		
	}
	
	public static double round(double value, int places) {
	    if (places < 0) throw new IllegalArgumentException();

	    long factor = (long) Math.pow(10, places);
	    value = value * factor;
	    long tmp = Math.round(value);
	    return (double) tmp / factor;
	}
	
	public List<Integer> twoOpt(List<Integer> perm, int i, int j){
        ArrayList<Integer> newPerm = new ArrayList<>(perm.subList(0, i));
        ArrayList<Integer> reversedPortion =  new ArrayList<>(perm.subList(i, j + 1));
        Collections.reverse(reversedPortion);
        newPerm.addAll(reversedPortion);
        newPerm.addAll(perm.subList(j + 1, perm.size()));

        return newPerm;
    }
	
	public double euc2d(double[] c1, double[] c2){
        return round(Math.sqrt(Math.pow(c1[0] - c2[0], 2.0) + Math.pow(c1[1] - c2[1], 2.0)), 0);
    }
	
	public double cost(List<Integer> permutation, double[][] cities){
        double distance = 0;
        for (int i = 0; i < permutation.size(); i++) {
            int c1 = i;
            int c2 = (i == permutation.size()-1)? permutation.get(0) : permutation.get(i+1);
            distance += euc2d(cities[c1], cities[c2]);
        }
        return round(distance, 4);
    }
	
	public List<Integer> stochasticTwoOpt(List<Integer> perm){
        Random rand = new Random();
		int c1 = rand.nextInt(perm.size());
        while (c1 == 0)
            c1 = rand.nextInt(perm.size());
        int c2 = rand.nextInt(perm.size());
        ArrayList<Integer> exclude = new ArrayList<>(Arrays.asList(c1, 0));
        exclude.add((c1==0) ? perm.size()-1 : c1-1);
        exclude.add((c1 == (perm.size()-1)) ? 0 : c1+1);
        while (exclude.contains(c2))
            c2 = rand.nextInt(perm.size());
        if (c2 < c1){
            int temp = c1;
            c1 = c2;
            c2 = temp;
        }
        return twoOpt(perm, c1, c2);
	}
	
	public Candidate localSearch(Candidate best, double[][] cities, int max_no_improv){
		int count = 0;
		while(count < max_no_improv) {
			Candidate candidate = new Candidate();
			candidate.vector = stochasticTwoOpt(best.vector);
			candidate.cost = cost(candidate.vector, cities);
			count = (candidate.cost < best.cost) ? 0:count+1;
			if(best.cost > candidate.cost){
				best = candidate;
			}
		}
		return best;
	}
	
	public Candidate contructRandomizedGreedySolution(double[][] cities, float alpha) {
		Random rand = new Random();
		Candidate cadidate = new Candidate();
		cadidate.vector.add(rand.nextInt(cities.length)); 
		
		List<Integer> allCities = new ArrayList<Integer>();
		for(int i = 0; i < cities.length; i++){
			allCities.add(i);
		}
		while(cadidate.vector.size() < cities.length) {
			List<Integer> candidates = new ArrayList<Integer>();
			for(int i=0; i <allCities.size(); i++){
				if(!cadidate.vector.contains(allCities.get(i))){
					candidates.add(allCities.get(i));
				}
			}
			List<Double> costs = new ArrayList<Double>();
			for(int i=0; i<candidates.size(); i++){
				double temp = euc2d(cities[cadidate.vector.size()-1], cities[candidates.get(i)]);
				costs.add(temp);
			}
			List<Integer> rcl = new ArrayList<Integer>();
			double min = Collections.min(costs);
			double max = Collections.max(costs);
			
			for(int i=0; i < costs.size(); i++){
				if(costs.get(i) <= (min + alpha*(max - min))){
					rcl.add(candidates.get(i));
				}
			}
			cadidate.vector.add(rcl.get(rand.nextInt(rcl.size())));
		}
		cadidate.cost = cost(cadidate.vector, cities);
		return cadidate;
	}
	
	public Candidate search(double[][] cities, int max_iter, int max_no_improve, float alpha) {
		Candidate best = new Candidate();
		best.cost = Double.MAX_VALUE;
		for(int i = 0; i < max_iter; i++){
			Candidate cadidate = contructRandomizedGreedySolution(cities, alpha);
			cadidate = localSearch(cadidate, cities, max_no_improve);
			if(best.cost > cadidate.cost){
				best = cadidate;
			}
			if (i % 1000 == 0) {
				double percent = round((double) i /max_iter * 100, 1);
				String progress = "\rProgress " + Double.toString(percent) + "%:\tCurrent best value:" + best.cost;
				System.out.print(progress);
			}
		}
		
		
		
		return best;
	}

	public double[][] parseString2Data(String input) {
		// TODO Auto-generated method stub
		return null;
	}
	
	public static double[][] berlin52 = {{565,575},{25,185},{345,750},{945,685},{845,655},
            {880,660},{25,230},{525,1000},{580,1175},{650,1130},{1605,620},
            {1220,580},{1465,200},{1530,5},{845,680},{725,370},{145,665},
            {415,635},{510,875},{560,365},{300,465},{520,585},{480,415},
            {835,625},{975,580},{1215,245},{1320,315},{1250,400},{660,180},
            {410,250},{420,555},{575,665},{1150,1160},{700,580},{685,595},
            {685,610},{770,610},{795,645},{720,635},{760,650},{475,960},
            {95,260},{875,920},{700,500},{555,815},{830,485},{1170,65},
            {830,610},{605,625},{595,360},{1340,725},{1740,245}};

	public static void main(String[] args) {
		// TODO Auto-generated method stub
//		int a = 5;
//		int b = a << 2;
//		int[] c = {1, 2};
//		int[] d = c;
//		Random rand = new Random();
//		int  n = rand.nextInt(50) + 1;
//		//System.out.println(round(200.23123, 0));
//		
//		int[][] f = {{1, 2}, {1, 2}};
//		int[][] f1 = f;
//		System.out.println(f1[1][1]);
		
		
		double[][] cities = berlin52;
		
		// Display data
		System.out.println("/***************************************************************/");
		System.out.println("The Greedy Randomized Adaptive Search - GRASP");
		System.out.println("The Berlin52 instance for GRASP");
		System.out.println("The optimal tour distance for Berlin52 instance is 7542 units");
		System.out.println("Author: Trunggm");
		System.out.println("/***************************************************************/");
		System.out.println("Data:");
		System.out.println("* Number of cities: " + cities.length);
		System.out.println("* List coordinates cities: ");
		for (int i = 0; i < cities.length; i+=4) {
			System.out.println(i + "/ [" + cities[i][0] + ", " + cities[i][1] + "]\t"
								+ (i+1) + "/ [" + cities[i+1][0] + ", " + cities[i+1][1] + "]\t"
								+ (i+2) + "/ [" + cities[i+2][0] + ", " + cities[i+2][1] + "]\t"
								+ (i+3) + "/ [" + cities[i+3][0] + ", " + cities[i+3][1] + "]");
		}
		
		
		int max_iter = 100000;
		int max_no_improve = 50;
		float greediness_factor = (float) 0.3;
		
		System.out.println("/***************************************************************/");
		System.out.println("Compute:");
		System.out.println("* Max iteration: \t\t" + max_iter);
		System.out.println("* Max number of improve: \t" + max_no_improve);
		System.out.println("* Greediness coefficient: \t" +  greediness_factor);
		
		Compute com = new Compute();
		
		Candidate best = com.search(cities, max_iter, max_no_improve, greediness_factor);
		System.out.println("");
		System.out.println("/***************************************************************/");
		String result = "";
		
		for(int i=0; i<best.vector.size();i++){
			result += best.vector.get(i) + ((i<best.vector.size()-1) ? "->" : "");
		}
		
		
		System.out.println("Final:");
		System.out.println("* Solution computed: " + best.cost + " units.");
		System.out.println("* Final tour path: ");
		System.out.println(result);
		
	}
}