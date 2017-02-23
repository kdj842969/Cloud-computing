import java.io.IOException;
import java.util.*;

import org.apache.hadoop.fs.Path;
import org.apache.hadoop.conf.*;
import org.apache.hadoop.io.*;
import org.apache.hadoop.mapred.*;
import org.apache.hadoop.util.*;

public class InvertedIndex {

    public static class Map extends MapReduceBase implements Mapper<LongWritable, Text, Text, Text> {
        JobConf conf;
        public void configure (JobConf job){
            this.conf = job;
        }
        public void map (LongWritable docId, Text value, OutputCollector<Text, Text> output, Reporter reporter) throws IOException {
            //retrieve # keywords from JobConf
            int argc = Integer.parseInt (conf.get ("argc"));
            //get the current file name
            FileSplit fileSplit = (FileSplit) reporter.getInputSplit();
            String filename = "" + fileSplit.getPath().getName();
            //read each line of file split
            String line = value.toString();
            //tokenize each word with a space
            StringTokenizer tokenizer = new StringTokenizer(line);
            //check if it is one of the given keywords
            HashMap<String, Integer> keywords = new HashMap<String, Integer>();
            int count = 0;
            for(int i=0; i<argc; i++){
                keywords.put(conf.get("keyword"+i), count);
            }
            String token = new String();
            int sum = 0;
            while (tokenizer.hasMoreTokens()) {
                token = tokenizer.nextToken();
                if(keywords.containsKey(token)){
                    output.collect(new Text(token), new Text(filename));
                }
            }
            //generate and pass pairs of a keyword and a document id to Reduce
        }
    }

    
    //receive one of the keywords and all document ids
    public static class Reduce extends MapReduceBase implements Reducer<Text, Text, Text, Text> {
        private Text result=new Text();
        private Text value = new Text();
        public void reduce(Text key, Iterator<Text> values, OutputCollector<Text, Text> output, Reporter reporter) throws IOException {
            HashMap m = new HashMap();
            int count = 0;
            while(values.hasNext()){
                String str = values.next().toString();
                if((m!=null)&&(m.get(str)!=null)){
                    count = (int)m.get(str);
                    m.put(str, ++count);
                } else {
		     m.put(str, 1);
                }
            }
            output.collect(key, new Text(m.toString()));
        }
    }
  
    public static void main(String[] args) throws Exception {
        //input format:
        //hadoop jar invertedindexes.jar InvertedIndexes input output keyword1 keyword2 ...
        JobConf conf = new JobConf(InvertedIndex.class);
        conf.setJobName("InvertedIndex");
	
        conf.setOutputKeyClass(Text.class);
	    conf.setOutputValueClass(Text.class);
	
	    conf.setMapperClass(Map.class);
	    conf.setCombinerClass(Reduce.class);
	    conf.setReducerClass(Reduce.class);
	
	    conf.setInputFormat(TextInputFormat.class);
	    conf.setOutputFormat(TextOutputFormat.class);
	
	    FileInputFormat.setInputPaths(conf, new Path(args[0]));
	    FileOutputFormat.setOutputPath(conf, new Path(args[1]));
	
        conf.set("argc", String.valueOf(args.length-2));
        for(int i=0; i<args.length-2;i++){
            conf.set("keyword"+i, args[i+2]);
        }
        
        JobClient.runJob(conf);
    }
}
