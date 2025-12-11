package egps.module.ncov.analysis.egps.module.million.tree.treebuilder.multithread;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

/** 
 *
 * @Description: 分配行数
 *
 *
 * @createdDate Sep 23, 2020
 * @lastModifiedDate Sep 23, 2020
 * @author "yudalang"
 */
public class RowLocationAssigner {
	
	/**
	 * 
	 * @Title: allocateRegion   
	 * @Description: 计算下三角形的平均分配位点。
	 *   
	 * @param: @param n
	 * @param: @param divideTimes
	 * @param: @return      
	 * @return: List<int[]>      
	 * @throws
	 */
	public List<int[]> allocateRegion(int n, int divideTimes) {
		List<int[]> ret = new ArrayList<>();
		
		if (divideTimes > 1) {
			
			int middleUpperTrianglePoint = middleUpperTrianglePoint(n);
			
			List<int[]> upperTri = allocateRegion(middleUpperTrianglePoint, --divideTimes);
			ret.addAll(upperTri);
			List<int[]> allocateDownTrapezoid = allocateDownTrapezoid(middleUpperTrianglePoint, n, --divideTimes);
			ret.addAll(allocateDownTrapezoid);
		}else if (divideTimes == 1) {
			int middleUpperTrianglePoint = middleUpperTrianglePoint(n);
			List<int[]> upperTriangle = 
					Arrays.asList(new int[] {0,middleUpperTrianglePoint},
					new int[] {middleUpperTrianglePoint,n});
			
			ret.addAll(upperTriangle);
		}else{
			ret.add(new int[] {0,n});
		}
		
		return ret;
	}
	
	public List<int[]> allocateDownTrapezoid(int n1,int n2, int divideTimes) {
		List<int[]> ret = new ArrayList<int[]>();
		
		if (divideTimes > 1) {
			int middleDownTrapezoidPoint = middleDownTrapezoidPoint(n1, n2);
			int divideTimes2 = --divideTimes;
			ret.addAll(allocateDownTrapezoid(n1, middleDownTrapezoidPoint, divideTimes2));
			ret.addAll(allocateDownTrapezoid(middleDownTrapezoidPoint, n2, divideTimes2));
		}else {
			int t = (int) (Math.sqrt( 0.5 * (n2 * n2 + n1 * n1 + n2 - n1) ));
			List<int[]> asList = Arrays.asList(new int[] {n1,t},new int[] {t,n2});
			ret.addAll(asList);
		}
		
		return ret;
	}
	
	/**
	 * 
	 * @Title: divideDownTrapezoid   
	 * @Description: divide the trapezoid
	 *   
	 * @param: @param n1
	 * @param: @param n2 should better than n1
	 * @param: @return      
	 * @return: List<int[]>      
	 */
	public List<int[]> divideDownTrapezoid(int n1,int n2) {
		int t = (int) (Math.sqrt( 0.5 * (n2 * n2 + n1 * n1 + n2 - n1) ));
		return Arrays.asList(new int[] {n1,t},new int[] {t,n2});
	}
	public int middleDownTrapezoidPoint(int n1,int n2) {
		int t = (int) (Math.sqrt( 0.5 * (n2 * n2 + n1 * n1 + n2 - n1) ));
		return t;
	}
	public List<int[]> divideUpperTriangle(int n) {
		int t = (int) (Math.sqrt(0.5 * n *(n - 1) + 0.25) + 0.5);
		return Arrays.asList(new int[] {0,t},new int[] {t,n});
	}
	public int middleUpperTrianglePoint(int n) {
		int t = (int) (Math.sqrt(0.5 * n *(n - 1) + 0.25) + 0.5);
		return t;
	}
	
	
	
	public static void main(String[] args) {
		RowLocationAssigner rowLocationAssigner = new RowLocationAssigner();
		List<int[]> allocateRegion = rowLocationAssigner.allocateRegionPrecise(8000,3);
		// should be 0,7 7 10
		System.out.println("Results: "+allocateRegion.size());
		for (int[] is : allocateRegion) {
			
			long a = is[0];
			long b = is[1];
			
			long totalNumber = (b - a + 1)*(a + b) /2;
			System.out.println(Arrays.toString(is) +"\t, Total number is :"+totalNumber);
		}
	}

	/**
	 * 
	 * @Title: allocateRegionPrecise   
	 * @Description: 不利用公式，比较精确地分割下三角！
	 *   注意，如果 numberOfChunks < 2或者
	 *   n < 101 就不会分割，返回一条线程的结果。
	 * @param: n : 总数
	 * @param: numberOfChunks : 目标区块地数量
	 * @return: List<int[]>      
	 */
	public List<int[]> allocateRegionPrecise(int n, int numberOfChunks) {
		LinkedList<int[]> ret = new LinkedList<int[]>();
		if (numberOfChunks < 2 || n < 101) {
			ret.add(new int[] {0,n});
			return ret;
		}
		
		long totalNumber = n;
		totalNumber = totalNumber * (totalNumber - 1) / 2;
		long interval = totalNumber / numberOfChunks;
		//System.out.println("total\t"+totalNumber + "\tinterval\t"+interval);
		
		boolean isFirst = true;
		
		long sum = 0;
		for (int i = 1; i <= n; i++) {
			sum += i;
			if (sum >= interval) {
				//System.out.println("达到interval, sum is "+sum + " i is: "+i);
				if (isFirst) {
					ret.add(new int[] {0,i});
					isFirst = false;
				}else {
					int prevSecond = ret.getLast()[1];
					ret.add(new int[] {prevSecond,i});
				}
				
				sum = 0;
			}
		}
		
		int prevSecond = ret.getLast()[1];
		
		if (prevSecond != n ) {
			ret.add(new int[] {prevSecond,n});
		}
		
		
		return ret;
	}
}
