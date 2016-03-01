import java.io.FileInputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;

/**
 * This is for COMP 3331/9331: Computer Networks & Applications Programming
 * Assignment: Routing Performance Analysis
 * 
 * @due_date: Friday 24 Oct 2014, 11:59 pm (Week 13)
 * @author Xi Zhang(z3472528), Xu Yan()
 *
 */
public class RoutingPerformance {
	private static ArrayList<Connection> requests = new ArrayList<Connection>(); // original
																					// connections
	private static List<Integer>[][] table;// routing table
	// correspondence relationship between router name and sequential number.
	private static Hashtable<String, Integer> routers = new Hashtable<String, Integer>();
	private static Scanner s;// topology.txt
	private static Scanner s1;// topology.txt
	private static Scanner s2;// workload.txt
	private static int num_of_connections = 0;
	private static int num_of_sucessful_connections = 0;
	private static int num_of_routers = 0;
	private static int num_of_packets = 0;
	private static int num_of_hops = 0;
	private static int packetLost = 0;
	private static int totalDelay = 0;
	private static int[][] hops; // adjacent relationship of routers on hops
	private static int[][] delay;// adjacent relationship of routers on delays
	private static int[][] load;// adjacent relationship of routers on capacity
	private static int[][] firm_load;// copy of firm_load
	private static float[][] occupation;// adjacent relationship of routers on
										// occupation

	/**
	 * Main Function
	 * 
	 * @param args
	 * @throws Exception
	 */
	@SuppressWarnings("unchecked")
	public static void main(String[] args) throws Exception {

		if (args.length != 5) {
			System.out.println("Required 5 arguments");
			return;
		}
		String NETWORK_SCHEME = args[0]; // CIRCUIT or PACKET
		String ROUTING_SCHEME = args[1]; // SHP or SDP or LLP
		String TOPOLOGY_FILE = args[2];
		String WORKLOAD_FILE = args[3];
		int PACKET_RATE = Integer.parseInt(args[4]);
		/*
		 * rewrite the comparator, so that it can match the requirement of
		 * collection sort on a particular variable.
		 */
		Comparator<Connection> c = new Comparator<Connection>() {
			@Override
			public int compare(Connection o1, Connection o2) {
				if (o1.start_time < o2.start_time) {
					return -1;
				}
				return 1;
			}
		};
		s = new Scanner(new FileInputStream(TOPOLOGY_FILE));
		s1 = new Scanner(new FileInputStream(TOPOLOGY_FILE));
		s2 = new Scanner(new FileInputStream(WORKLOAD_FILE));
		/* first scan topolofy.txt to store the router names in a hashtable. */
		while (s.hasNextLine()) {
			String line = s.nextLine();
			String[] temp = line.split(" ");
			routers.put(temp[0], 1);
			routers.put(temp[1], 1);
		}
		Set<String> keys = routers.keySet();
		/*
		 * we numbered the routers stored in the hashtable by a random sequence
		 * but with an group of consecutive numbers starting from 0.(notice:
		 * this matches every kind of situations except when more than one
		 * routers have the same router name. The correspondence between router
		 * names and numbers are 1-1 but router names are not in alphabetical
		 * order. This will leads to different results from those which uses
		 * other orders of storing data when SHP schemes are used.)
		 */
		for (String key : keys) {
			routers.put(key, num_of_routers++);
		}
		hops = new int[num_of_routers][num_of_routers];
		delay = new int[num_of_routers][num_of_routers];
		load = new int[num_of_routers][num_of_routers];
		firm_load = new int[num_of_routers][num_of_routers];
		occupation = new float[num_of_routers][num_of_routers];
		for (int i = 0; i < num_of_routers; i++) {
			for (int j = 0; j < num_of_routers; j++) {
				hops[i][j] = 99999;
				delay[i][j] = 99999;
				load[i][j] = 0;
				firm_load[i][j] = 99999;
				occupation[i][j] = 0;
			}
		}
		/*
		 * second scan topology.txt to store the data in to different
		 * containers(hops,delay,load,firm_load)
		 */
		while (s1.hasNextLine()) {
			String line = s1.nextLine();
			String[] temp = line.split(" ");
			hops[routers.get(temp[0])][routers.get(temp[1])] = 1;
			delay[routers.get(temp[0])][routers.get(temp[1])] = Integer
					.parseInt(temp[2]);
			firm_load[routers.get(temp[0])][routers.get(temp[1])] = Integer
					.parseInt(temp[2]);
			load[routers.get(temp[0])][routers.get(temp[1])] = Integer
					.parseInt(temp[3]);
			hops[routers.get(temp[1])][routers.get(temp[0])] = 1;
			delay[routers.get(temp[1])][routers.get(temp[0])] = Integer
					.parseInt(temp[2]);
			firm_load[routers.get(temp[1])][routers.get(temp[0])] = Integer
					.parseInt(temp[2]);
			load[routers.get(temp[1])][routers.get(temp[0])] = Integer
					.parseInt(temp[3]);
		}
		/*
		 * scan workload into a new class connection then connect them with a
		 * ArrayList then sort them against starting time.
		 */
		while (s2.hasNextLine()) {
			String line = s2.nextLine();
			String[] temp = line.split(" ");
			Connection e = new Connection(Double.valueOf(temp[0]),
					Double.valueOf(temp[3]) + Double.valueOf(temp[0]),
					routers.get(temp[1]), routers.get(temp[2]));
			requests.add(e);
		}
		Collections.sort(requests, c); // sorting
		table = new List[num_of_routers][num_of_routers];
		/*
		 * Using Dijkstra Algorithm to calculate the whole routing table
		 * according to different routing schemes. Since we using SDP scheme to
		 * calculate first routing table when LLP is chosen, SDP and LLP schemes
		 * share the same code below.
		 */
		if (ROUTING_SCHEME.equals("SHP")) {
			for (int start = 0; start < num_of_routers; start++) {
				int[][] tempHops = new int[num_of_routers][num_of_routers];
				for (int i = 0; i < num_of_routers; i++) {
					for (int j = 0; j < num_of_routers; j++) {
						tempHops[i][j] = hops[i][j];
					}
				}
				Dijkstra dij = new Dijkstra(tempHops, start);
				for (int end = 0; end < num_of_routers; end++) {
					if (start == end) {
						continue;
					}
					Map<Integer, List<Integer>> pathList = dij.getPathListMap();
					table[start][end] = pathList.get(end);
				}
			}
		} else { // SDP or LLP(first time)
			for (int start = 0; start < num_of_routers; start++) {
				int[][] tempDelay = new int[num_of_routers][num_of_routers];
				for (int i = 0; i < num_of_routers; i++) {
					for (int j = 0; j < num_of_routers; j++) {
						tempDelay[i][j] = delay[i][j];
					}
				}
				Dijkstra dij = new Dijkstra(tempDelay, start);
				for (int end = 0; end < num_of_routers; end++) {
					if (start == end) {
						continue;
					}
					Map<Integer, List<Integer>> pathList = dij.getPathListMap();
					table[start][end] = pathList.get(end);
				}
			}
		}
		/* processing packets transmission here */
		if (NETWORK_SCHEME.equals("CIRCUIT")) {
			circuitNetwork(ROUTING_SCHEME, PACKET_RATE);
		} else { // NETWORK_SCHEME.equals("PARCKET")
			packetNetwork(ROUTING_SCHEME, PACKET_RATE);
		}
		/************************** Outputs ***********************************/
		float blockRate = (float) packetLost / num_of_packets;
		System.out.printf("total number of virtual circuit requests: "
				+ "%d\n\n", num_of_connections);
		System.out.printf("total number of packets: " + "%d\n\n",
				num_of_packets);
		System.out.printf("number of successfully routed packets: " + "%d\n\n",
				num_of_packets - packetLost);
		System.out.printf("percentage of successfully routed packets: "
				+ "%.2f\n\n", 100 - blockRate * 100);
		System.out.printf("number of blocked packets: " + "%d\n\n", packetLost);
		System.out.printf("percentage of blocked packets: " + "%.2f\n\n",
				blockRate * 100);
		System.out.printf("average number of hops per circuit: " + "%.2f\n\n",
				(float) num_of_hops / num_of_sucessful_connections);
		System.out.printf("average cumulative propagation delay per circuit: "
				+ "%.2f\n", (float) totalDelay / num_of_sucessful_connections);
	}

	/**
	 * no pre-processes on CIRCUIT network scheme
	 * 
	 * @param scheme
	 * @param n
	 */
	private static void circuitNetwork(String scheme, int n) {
		process(requests, scheme, n);
	}

	/**
	 * pre-process deviding different packets into different connections by
	 * giving different starting time and ending time to them, then sorting them
	 * against new starting time
	 * 
	 * @param scheme
	 * @param n
	 */
	private static void packetNetwork(String scheme, int n) {
		ArrayList<Connection> newReq = new ArrayList<Connection>();
		double duration = (double) 1 / n;
		Iterator<Connection> it = requests.iterator();
		while (it.hasNext()) {
			Connection unit = it.next();
			int times = (int) (int) (n * (unit.end_time - unit.start_time) + 1);
			for (int l = 0; l < times; l++) {
				double real_start = unit.start_time + l * duration;
				double tmp_real_end = unit.start_time + (l + 1) * duration;
				double real_end = (tmp_real_end > unit.end_time) ? unit.end_time
						: tmp_real_end;
				Connection e = new Connection(real_start, real_end,
						unit.source, unit.destination);
				newReq.add(e);
			}
		}
		Comparator<Connection> c = new Comparator<Connection>() {
			@Override
			public int compare(Connection o1, Connection o2) {
				if (o1.start_time < o2.start_time) {
					return -1;
				}
				return 1;
			}
		};
		Collections.sort(newReq, c);
		process(newReq, scheme, 1);
	}

	/**
	 * After pre-processing (no processing under CIRCUIT scheme), real
	 * transmission tasks are carried out in the method below. There is no
	 * difference in this method between CIRCUIT and PACKET schemes. but there
	 * are differences between SDP/SHP and LLP. SDP/SHP only needs to use the
	 * original routing table generated first time by their certain adjacent
	 * relationship and costs in between. As LLP is a dynamic routing algorithm,
	 * the routing table changes each time a connection is setup.
	 * 
	 * @param connection
	 * @param scheme
	 * @param rate
	 */
	private static void process(ArrayList<Connection> connection,
			String scheme, int rate) {
		ArrayList<BusyPaths> busyPath = new ArrayList<BusyPaths>();
		try {
			int initial = 0; // a flag
			for (Connection con : connection) {
				num_of_connections++;
				num_of_packets += (int) (rate * (con.end_time - con.start_time) + 1);
				List<Integer> route; // route which stores the current most
										// optimized route
				if (scheme.equals("SHP") || scheme.equals("SDP")) {
					route = table[con.getSource()][con.getDestination()];
				} else {
					/*
					 * using routing table generated by SDP scheme for the first
					 * time.
					 */
					if (initial == 0) {
						route = table[con.getSource()][con.getDestination()];
						initial++;
					} else {
						int start = con.getSource();
						float[][] tempOccup = new float[num_of_routers][num_of_routers];
						for (int i = 0; i < num_of_routers; i++) {
							for (int j = 0; j < num_of_routers; j++) {
								occupation[i][j] = 1 - (float) load[i][j]
										/ firm_load[i][j];
								tempOccup[i][j] = occupation[i][j];
							}
						}
						LLPDijkstra dij = new LLPDijkstra(tempOccup, start);
						for (int end = 0; end < num_of_routers; end++) {
							if (start == end) {
								continue;
							}
							Map<Integer, List<Integer>> pathList = dij
									.getPathListMap();
							table[start][end] = pathList.get(end);
						}
						route = table[con.getSource()][con.getDestination()];
					}
				}
				int pathLength = route.size();
				Iterator<Integer> r = route.iterator();
				int[] mark = new int[pathLength]; // using this array to store
													// the mark of the route
				int k = 0;
				while (r.hasNext()) {
					mark[k++] = r.next();
				}
				release(busyPath, con); // release the busy path
				int lost = 0; // a flag indicating when packet lost happens (1)
				/*
				 * if there is one part of the whole route whose capacity is 0.
				 * the whole connection request will be discarded.
				 */
				for (int i = 0, j = 1; i < pathLength - 1 && j < pathLength; i++, j++) {
					if (load[mark[i]][mark[j]] == 0
							|| load[mark[j]][mark[i]] == 0) {
						lost = 1;
						break;
					}
				}
				/*
				 * After confirming that there is no 0 capacity along side the
				 * whole path, the current path will be marked as "busy" and
				 * store every single section of the whole path into busy path
				 * ArrayList.
				 */
				if (lost == 0) {
					num_of_hops += (pathLength - 1);
					num_of_sucessful_connections++;
					for (int i = 0, j = 1; i < pathLength - 1 && j < pathLength; i++, j++) {
						totalDelay += delay[mark[i]][mark[j]];
						/*
						 * In a two-way network, capacity is shared by each end.
						 */
						load[mark[i]][mark[j]]--;
						load[mark[j]][mark[i]]--;
						BusyPaths temp1 = new BusyPaths(con.start_time,
								mark[i], mark[j], con.end_time);
						BusyPaths temp2 = new BusyPaths(con.start_time,
								mark[j], mark[i], con.end_time);
						busyPath.add(temp1);
						busyPath.add(temp2);
					}
				} else {
					packetLost += (int) (rate * (con.end_time - con.start_time) + 1);
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	/**
	 * if the starting time of current connection is later than the ending time
	 * of busy path in the busyPath list, these busy path will be set idle and
	 * capacity will be released by incrementing 1 per release.
	 * 
	 * @param busy
	 * @param c
	 */
	private static void release(ArrayList<BusyPaths> busy, Connection c) {
		for (int i = 0; i < busy.size(); i++) {
			try {
				BusyPaths single = busy.get(i);
				if (c.start_time >= single.getFinish()) {
					load[single.getSource()][single.getDestination()]++;
					i--;
					busy.remove(single);
				}
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}
}

class Dijkstra {
	private int[][] graph;
	private int start;
	private int dimention;
	private Map<Integer, List<Integer>> pathListMap = new HashMap<Integer, List<Integer>>();

	/**
	 * Constructor of Dijkstra class
	 * 
	 * @param graph
	 * @param start
	 */
	public Dijkstra(int[][] graph, int start) {
		this.graph = graph;
		this.start = start;
		this.dimention = graph.length;
		calculate();
	}

	/**
	 * calculate the shortest path according to a specific matrix
	 * 
	 */
	private void calculate() {

		for (int end = 0; end < dimention; end++) {
			List<Integer> pathList = new ArrayList<Integer>();
			pathList.add(start);
			pathList.add(end);
			pathListMap.put(end, pathList);
		}
		int[] visited = new int[dimention]; // store visited vertices
		visited[start] = 1;
		/* there are n - 1 vertices which needs traversal. */
		for (int count = 1; count < dimention; count++) {
			int bridge = -1;
			int cost = Integer.MAX_VALUE;
			for (int end = 0; end < dimention; end++) {
				if (visited[end] == 0 && graph[start][end] < cost) {
					cost = graph[start][end];
					bridge = end;
				}
			}
			visited[bridge] = 1;
			/*
			 * if distance between start and end stored in graph now is bigger
			 * than the plus of distances from the start to the bridge and from
			 * the bridge to the end, then replace the old one with the new cost
			 * which is less than before.
			 */
			for (int end = 0; end < dimention; end++) {
				if (visited[end] == 0
						&& graph[start][bridge] + graph[bridge][end] < graph[start][end]) {
					graph[start][end] = graph[start][bridge]
							+ graph[bridge][end];
					List<Integer> pathList = pathListMap.get(end);
					List<Integer> bridgePathList = pathListMap.get(bridge);
					pathList.clear();
					pathList.addAll(bridgePathList);
					pathList.add(end);
				}
			}
		}
	}

	/**
	 * @return
	 */
	public Map<Integer, List<Integer>> getPathListMap() {
		return pathListMap;
	}
}

class LLPDijkstra {
	private float[][] graph;
	private int start;
	private int dimention;
	private Map<Integer, List<Integer>> pathListMap = new HashMap<Integer, List<Integer>>();

	/**
	 * @param graph
	 * @param start
	 */
	public LLPDijkstra(float[][] graph, int start) {
		this.graph = graph;
		this.start = start;
		this.dimention = graph.length;
		calculate();
	}

	/**
	 * 
	 */
	private void calculate() {

		for (int end = 0; end < dimention; end++) {
			List<Integer> pathList = new ArrayList<Integer>();
			pathList.add(start);
			pathList.add(end);
			pathListMap.put(end, pathList);
		}
		int[] visited = new int[dimention];
		visited[start] = 1;
		for (int count = 1; count < dimention; count++) {
			int bridge = -1;
			float cost = Integer.MAX_VALUE;
			for (int i = 0; i < dimention; i++) {
				if (visited[i] == 0 && graph[start][i] < cost) {
					cost = graph[start][i];
					bridge = i;
				}
			}
			visited[bridge] = 1;
			/*
			 * if load between start and end stored in graph now is bigger
			 * than the maximum loads of each route from the start to the bridge and from
			 * the bridge to the end, then replace the old load with the new load
			 * which is less than before.
			 */
			for (int i = 0; i < dimention; i++) {
				if (visited[i] == 0
						&& Math.max(graph[start][bridge], graph[bridge][i]) < graph[start][i]) {
					graph[start][i] = Math.max(graph[start][bridge],
							graph[bridge][i]);
					List<Integer> pathList = pathListMap.get(i);
					List<Integer> bridgePathList = pathListMap.get(bridge);
					pathList.clear();
					pathList.addAll(bridgePathList);
					pathList.add(i);
				}
			}
		}
	}

	public float[][] getGraph() {
		return graph;
	}

	public int getStart() {
		return start;
	}

	public int getDimention() {
		return dimention;
	}

	public Map<Integer, List<Integer>> getPathListMap() {
		return pathListMap;
	}

}

class BusyPaths {

	private int source;
	private int destination;
	private double start;
	private double finish;

	public BusyPaths(double start, int src, int dest, double finish) {
		this.source = src;
		this.destination = dest;
		this.start = start;
		this.finish = finish;
	}

	public int getSource() {
		return source;
	}

	public int getDestination() {
		return destination;
	}

	public double getStart() {
		return start;
	}

	public double getFinish() {
		return finish;
	}

}

class Connection {

	double start_time;
	double end_time;
	int source;
	int destination;

	public Connection(double s, double d, int sou, int des) {
		this.start_time = s;
		this.end_time = d;
		this.source = sou;
		this.destination = des;
	}

	public double getStart_time() {
		return start_time;
	}

	public double getEnd_time() {
		return end_time;
	}

	public int getSource() {
		return source;
	}

	public int getDestination() {
		return destination;
	}

}
