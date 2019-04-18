use std::io::{BufRead, BufReader};
use std::fs::File;
use std::collections::HashMap;
use std::process;
use core::cmp;

use std::io::prelude::*;
use std::fs::OpenOptions;
use std::thread;
use std::fs;
use std::sync::{mpsc, Arc};
use std::env;

//const CHROMOSOMES: &'static [&'static str] = &["chr1","chr2","chr3","chr4","chr5","chr6","chr7",
//"chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
//"chr20","chr21","chr22", "chrX", "chrY"];

const CHROMOSOMES: &'static [&'static str] = &["chr1","chrY"];



fn main() {

    let args : Vec<String> = env::args().collect();
    let readpath = args.get(1).unwrap(); //"/home/menzel/Desktop/THM/promotion/projekte/nopeak/rust_profile/test.bed";

    let genomepath  = args.get(2).unwrap(); // "/home/menzel/hg19.fa";
    let q: i32 = args.get(3).unwrap().parse().unwrap(); //8;
    let radius : i32 =  args.get(4).unwrap().parse().unwrap(); //500;

    create_profile(readpath.to_string(), genomepath.to_string(), q, radius);
}

fn create_profile(readpath: String, genomepath: String, q : i32, radius : i32){


    println!("v {}", 5);
    // read data

    let mut tmp : (HashMap<String, Vec<i32>>, HashMap<String, Vec<i32>>);
    println!("{} from {}", "read reads", readpath);

    tmp = read_reads(&readpath);

    let reads_p = tmp.0;
    let reads_n = tmp.1;

    println!("{}", "read genome");
    let genome = Arc::new(read_chr(&genomepath));
    let g1 = Arc::clone(&genome);
    let g2 = Arc::clone(&genome);

    /////////////////////////////
    // multithread
    /////////////////////////////

    // create profiles for n and p

    let (s1, r1) = mpsc::channel();
    let (s2, r2) = mpsc::channel();

    println!("{}", "build profiles");
    let mut prof_n_folded: HashMap<String, Vec<i32>> = HashMap::new();
    let mut prof_p_folded: HashMap<String, Vec<i32>> = HashMap::new();


    let n_handle= thread::spawn( move || {
        s1.send(fold_profile(profile(reads_n, &g1, q, radius))).expect("something wrong with the profile");
    });

    let p_handle= thread::spawn( move || {
        s2.send(fold_profile(profile(reads_p, &g2, q, radius))).expect("something wrong with the profile");
    });

    n_handle.join().expect("something wrong with the thread");
    p_handle.join().expect("something wrong with the thread");

    /////////////////////////////
    // multithread
    /////////////////////////////


    let t1 = r1.recv().unwrap();

    for k in t1.keys(){
        prof_n_folded.insert(k.clone(), t1.get(k).unwrap().clone());
    }

    let t2 = r2.recv().unwrap();

    for k in t2.keys(){
        prof_p_folded.insert(k.clone(), t2.get(k).unwrap().clone());
    }

    // merge profile
    println!("{}", "merge profiles");
    let result = merge_fw_and_bw(prof_p_folded, prof_n_folded);

    //write profile to file
    println!("{}", "write profiles");
    let filepath = "/tmp/result";
    let readcount = 0;
    write_to_file(filepath, result, genomepath.to_string(), readpath.to_string(), radius, q, readcount);

}

fn merge_fw_and_bw(prof_p_folded: HashMap<String, Vec<i32>>, mut prof_n_folded: HashMap<String, Vec<i32>>) -> HashMap<String, Vec<i32>>{
    let mut updated : HashMap<String, Vec<i32>> = HashMap::new();


    for qmer in prof_p_folded.keys() {
        let qmer_rc: String;

        if prof_n_folded.contains_key(qmer) {
            qmer_rc = qmer.clone();
        } else if prof_n_folded.contains_key(&reverse_complement(qmer)) {
            qmer_rc = reverse_complement(qmer);
        } else {
            updated.insert(qmer.to_string(), prof_p_folded.get(qmer).unwrap().to_vec());
            continue;
        }

        let p1 = prof_p_folded.get(qmer).unwrap();
        let p2 = prof_n_folded.get(&qmer_rc).unwrap();

        let mut new: Vec<i32> = Vec::new();

        for v in 0..p1.len() {
            new.push(p1[v] + p2[p2.len() - 1 - v]);
        }
        new.reverse();

        prof_n_folded.remove(qmer);

        updated.insert(qmer.to_string(), new);
    }

    for qmer in prof_n_folded.keys(){
        let mut tmp = prof_n_folded.get(qmer).unwrap().to_vec();
        tmp.reverse();
        updated.insert(qmer.to_string(), tmp);
    }

    updated
}


fn write_to_file(pathstr: &str, result : HashMap<String, Vec<i32>>, genomepath : String, readpath : String, radius : i32, q : i32, readcount :i32){

    if fs::metadata(pathstr).is_ok() {
        fs::remove_file(pathstr).expect("could not delete old file");
    }
    let mut file = OpenOptions::new().create(true).write(true).append(true).open(pathstr).unwrap();

    let comment = format!("# genome {}\n# reads: {}\n# radius: {}\n# q: {}\n# reads used: {}", genomepath, readpath, radius, q, readcount);
    if let Err(e) = writeln!(file, r"{}", comment) {
        eprintln!("Couldn't write to file: {}", e);
    }

    let header = ["qmer\t".to_string(), (0..(radius * 2 + 2)).map(|x| x.to_string()).collect::<Vec<String>>().join("\t")].concat();
    if let Err(e) = writeln!(file, r"{}", header) {
        eprintln!("Couldn't write to file: {}", e);
    }

    // sort by kmer before writing
    //let mut keys : Vec<String> = result.clone().into_iter().map(|(k, _v)| k.to_string()).collect();
    //keys.sort();

    for qmer in result.keys(){

        let counts : Vec<String> = result.get(qmer).unwrap().iter().map(|x| x.to_string()).collect();
        let out = [qmer.to_uppercase(), "\t".to_string(), counts.join("\t").to_string()].concat();

        if let Err(e) = writeln!(file, r"{}", out) {
            eprintln!("Couldn't write to file: {}", e);
        }
    }
}


fn fold_profile(profile : HashMap<String, Vec<i32>>) -> HashMap<String, Vec<i32>>{

    let mut updated : HashMap<String, Vec<i32>> = HashMap::new();

    for qmer in profile.keys(){
        let qmer_rc = reverse_complement(&qmer);

        if qmer > &qmer_rc && profile.contains_key(&qmer_rc){
            continue;
        }

        if &qmer_rc != qmer && profile.contains_key(&qmer_rc){

            let p1 = profile.get(qmer).unwrap();
            let p2 = profile.get(&qmer_rc).unwrap();

            let mut new: Vec<i32> = Vec::new();

            for v in 0..p1.len(){
                new.push(p1[v] + p2[v]);
            }

            updated.insert(qmer.to_string(), new);


        } else {

            updated.insert(qmer.to_string(), profile.get(qmer).unwrap().to_vec());
        }
    }

    let result = updated;

    result
}


fn reverse_complement(qmer : &str) -> String {

    let tmp : String = qmer.chars().rev().collect();

    tmp.replace("a", "x").replace("t", "a").replace("x","t").replace("g","x").replace("c","g").replace("x","c")
}

fn profile(reads: HashMap<String, Vec<i32>>, genome: &Arc<HashMap<String, String>>, q : i32, radius : i32) -> HashMap<String, Vec<i32>>{

    let mut profile : HashMap<String, Vec<i32>> = HashMap::new();

    for chr in CHROMOSOMES{
        println!("currently on {}", chr);
        let mut tmp = reads.get(*chr).unwrap();
        let mut stop = 0;

        if genome.contains_key(*chr) {
            let all = genome.get(*chr).unwrap();
            let chrl = all.len() as i32;

            let mut lastq: Vec<String> =  Vec::new();
            let mut lastd: Vec<i32> =  Vec::new();
            let mut has_neigbour = false;
            let mut diff = 0;

            for (it,read) in tmp.iter().enumerate() {

                stop += 1;

                if stop > 1{
                    process::exit(0x0);
                }

                let mut left = cmp::max(0, read - radius - q / 2);
                let mut right = cmp::min(read + radius + q / 2, chrl) - q;
                let mut i = 0;

                if has_neigbour{
                    has_neigbour = false;
                    for (s,qmer) in lastq.iter().enumerate(){
                        if lastd[s] > 0 && lastd[s] < radius - diff {
                            let mut vec: &mut Vec<i32> = profile.entry(qmer.parse().unwrap()).or_default();
                            let pos = (diff + lastd[s]) as usize;
                            //println!("diff {} last {} pos {} {:?} {:?}", diff, lastd[s], pos ,lastd, lastq);

                            vec[pos] = vec[pos] + 1;

                            println!("Setting +1 at {} for {} ",pos , qmer);
                        }
                    }

                    left = right - diff;
                    //println!("right {} new left {}", right , left );
                }

                if it + 1 < tmp.len() {
                    diff = tmp.get(it + 1).unwrap() - read;
                } else {
                    diff = 100;
                }

                if  diff < 10 { //if there is a next read and this one is less  than 10 pos away

                    if lastd.len() > 1 {
                        for v in 0..lastd.len() {
                            lastd[v] = lastd[v] - diff;
                        }
                    }


                    has_neigbour = true;
                    //println!("read {} neigbour {} diff {} ", read, tmp.get(it+1).unwrap(), diff);

                } else {
                    let mut lastq: Vec<String> =  Vec::new();
                    let mut lastd: Vec<i32> =  Vec::new();
                    left = cmp::max(0, read - radius - q / 2);
                }


                for p in left..right + 1{
                    let qmer = &all[p as usize..(p + (q as i32)) as usize];

                    if qmer.contains("n") {
                        i += 1;
                        continue;
                    }

                    let vec = profile.entry(qmer.to_string()).or_insert(vec![0; (radius * 2 + 1) as usize]);

                    vec[i] = vec[i] + 1;
                    i += 1;
                    println!("Setting +1 at {} for {} ",i ,qmer);

                    if has_neigbour && left > diff {
                        lastq.push(qmer.parse().unwrap());
                        lastd.push(i as i32);
                    }
                }


                println!("original:");

                right = cmp::min(read + radius + q / 2, chrl) - q;
                left = cmp::max(0, read - radius - q / 2);


                for p in left..right + 1{
                    let qmer = &all[p as usize..(p + (q as i32)) as usize];

                    if qmer.contains("n") {
                        i += 1;
                        continue;
                    }

                    println!("Setting +1 at {} for {} ",i ,qmer);
                }

            }
        }
    }


    profile
}


fn read_chr(genomepath: &String) -> HashMap<String, String>{
    let file = File::open(genomepath).unwrap();
    let mut lines : Vec<String> = Vec::new();
    let mut genome : HashMap<String, String> = HashMap::new();
    let mut lastchr: String = "".parse().unwrap();

    for line in BufReader::new(file).lines() {
        let l = line.unwrap();

        if l.starts_with(">chr"){

            if lines.len() > 0 {
                genome.insert(lastchr.clone()[1..].trim().to_string(), lines.concat().to_ascii_lowercase());
                lines.clear();
            }
            lastchr = l.clone();

        } else {
            lines.push(l);
        }
    }

    genome.insert(lastchr.clone()[1..].trim().to_string(), lines.concat().to_ascii_lowercase());

    genome
}


fn read_reads(path: &String) -> (HashMap<String, Vec<i32>>, HashMap<String, Vec<i32>>){

    let file = File::open(path).unwrap();
    let mut readstarts_n: HashMap<String, Vec<i32>> = HashMap::new();
    let mut readstarts_p: HashMap<String, Vec<i32>> = HashMap::new();

    let mut tmp_n : Vec<i32> = Vec::new();
    let mut tmp_p : Vec<i32> = Vec::new();
    let mut oldchr : String = "".parse().unwrap();
    let mut readcount = 0;

    for line in BufReader::new(file).lines() {
        let l = line.unwrap();
        let parts : Vec<&str> = l.split("\t").collect();

        if oldchr != parts[0]{

                if oldchr.len() > 1 {
                    tmp_n.sort();
                    tmp_n.dedup();
                    tmp_p.sort();
                    tmp_p.dedup();
                    readstarts_p.insert(oldchr.clone(), tmp_p.clone());
                    readstarts_n.insert(oldchr.clone(), tmp_n.clone());
                }

                oldchr = parts[0].to_string();

                tmp_n.clear();
                tmp_p.clear();
            }

        if parts[5] == "+" {
            tmp_p.push(parts[1].parse::<i32>().unwrap().clone())
        } else {
            tmp_n.push((parts[2].parse::<i32>().unwrap() - 1).clone())
        }
        readcount += 1;

    }

    tmp_n.sort();
    tmp_p.sort();

    readstarts_p.insert(oldchr.clone(), tmp_p.clone());
    readstarts_n.insert(oldchr.clone(), tmp_n.clone());

    println!("read {} reads", readcount);

    (readstarts_p, readstarts_n)
}

