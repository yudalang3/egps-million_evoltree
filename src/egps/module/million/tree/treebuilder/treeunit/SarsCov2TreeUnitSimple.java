package egps.module.ncov.analysis.egps.module.million.tree.treebuilder.treeunit;

import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Objects;

public class SarsCov2TreeUnitSimple implements TreeNodeUnit {

	final double dayAfterEarliestDate;

	final String country;

	Map<String, Double> country2weightMap;
	
	/**
	 * 国家数量阈值，可根据内存占用量适当调整
	 */
	int countryCountThreshold=5;

	public SarsCov2TreeUnitSimple(double dayAfterEarliestDate, String country, Map<String, Double> country2weightMap) {
		this.dayAfterEarliestDate = dayAfterEarliestDate;
		this.country = country;
		this.country2weightMap = country2weightMap;
	}

	public Map<String, Double> getCountry2weightMap() {
		return country2weightMap;
	}

	@Override
	public String toString() {
		return "day AfterEarliestDate is :" + dayAfterEarliestDate + " Country: " + country;
	}

	@Override
	public TreeNodeUnit produceCluster(TreeNodeUnit b) {

		SarsCov2TreeUnitSimple bCov2TreeUnit = (SarsCov2TreeUnitSimple) b;

		double tempDayAfterE = 0.5 * (this.dayAfterEarliestDate + bCov2TreeUnit.dayAfterEarliestDate);

		String country = null;

		// 新的用于给合并节点赋值的map
		Map<String, Double> tempMap = new HashMap<String, Double>();

		if (Objects.nonNull(this.country) && Objects.nonNull(bCov2TreeUnit.country)) {
			// 两个国家都不是空时，通过赋值map，然后取权重最大的作为国家

			// 先将当一个子节点的国家map按权重0.5,赋值给tempMap
			for (Entry<String, Double> entry : this.country2weightMap.entrySet()) {
				tempMap.put(entry.getKey(), entry.getValue() * 0.5);

			}

			// 将第二个国家的map按权重0.5加上去
			for (Entry<String, Double> entry : bCov2TreeUnit.country2weightMap.entrySet()) {
				String key = entry.getKey();
				Double value;
				if (tempMap.containsKey(key)) {
					value = tempMap.get(key) + entry.getValue() * 0.5;
				} else {
					value = entry.getValue() * 0.5;
				}
				tempMap.put(key, value);
			}

			// 如果国家数量超过一定阈值，现在设置的5，则认map留权重最高的5个;否则国家赋值权重最大的
			if (tempMap.size() > countryCountThreshold) {
				Map<String, Double> sortMap = sortByValue(tempMap);
				tempMap.clear();
				int count=0;
				for (Entry<String, Double> entry : sortMap.entrySet()) {
				      if(count<countryCountThreshold) {
				    	  tempMap.put(entry.getKey(), entry.getValue());
				    	  count++;
				      }
					
				}
				
			} else {
				String tempCountry = null;
				double tempWeight = 0;
				for (Entry<String, Double> entry : tempMap.entrySet()) {
					if (tempWeight < entry.getValue()) {
						tempWeight = entry.getValue();
						tempCountry = entry.getKey();
					}

				}
				country = tempCountry;
			}
			// =========简单判断的旧代码=============
//			if (this.country.equals(bCov2TreeUnit.country)) {
//				country = this.country;
//			}
		}

		SarsCov2TreeUnitSimple sarsCov2TreeUnit = new SarsCov2TreeUnitSimple(tempDayAfterE, country, tempMap);
		return sarsCov2TreeUnit;
	}

	/**
	 * 
	 * 让map按Value倒序排序
	 * 
	 * @title sortByValue
	 * @createdDate 2021-02-27 12:24
	 * @lastModifiedDate 2021-02-27 12:24
	 * @author yjn
	 * @since 1.7
	 * 
	 * @param <K>
	 * @param <V>
	 * @param map
	 * @return
	 * @return Map<K,V>
	 */
	public static <K, V extends Comparable<? super V>> Map<K, V> sortByValue(Map<K, V> map) {
		Map<K, V> result = new LinkedHashMap<>();

		map.entrySet().stream().sorted(Map.Entry.<K, V>comparingByValue().reversed())
				.forEachOrdered(e -> result.put(e.getKey(), e.getValue()));
		return result;
	}

	
	/**
	 * 
	 * 为了测试一下sort方法，后面可以移到eGPS的util里
	 *  
	 * @title main
	 * @createdDate 2021-02-27 12:31
	 * @lastModifiedDate 2021-02-27 12:31
	 * @author yjn
	 * @since 1.7
	 *   
	 * @param args
	 * @return void
	 */
	public static void main(String[] args) {
		Map<String, Double> map = new HashMap<>();
		map.put("one", 0.08);
		map.put("two", 0.1);
		map.put("three", 0.2);
		map.put("four", 0.91);
		for (Entry<String, Double> entry : map.entrySet()) {
			System.out.println(entry.getKey() + " " + entry.getValue());
		}

		Map<String, Double> sortMap = SarsCov2TreeUnitSimple.sortByValue(map);

		for (Entry<String, Double> entry : sortMap.entrySet()) {
			System.out.println(entry.getKey() + " " + entry.getValue());
		}

	}

}
