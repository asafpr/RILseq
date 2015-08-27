from subprocess import Popen,PIPE
from Queue import Queue
import threading
import sys
class Runner(threading.Thread):

    def __init__(self, task_queue, result_queue,RNAup):
        threading.Thread.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue
        self.p = Popen(RNAup, shell=True, stdin=PIPE,
                       stdout=PIPE, stderr=PIPE, close_fds=True)
        self.RNAup=RNAup
        
    def run(self):
        counter=0
        while True:
            next_task=self.task_queue.get()
            #next_task is an array: input sequence, srna, txid, target, the last three are used to 
            #in the return array
            if next_task is None:
                # Poison pill means shutdown
                try:
                    self.p.stdin.write('@\n')
                except:
                    pass
                self.task_queue.task_done()
                break
            inputseq = next_task[0]
            get_fold = next_task[1] == 'get_fold'
            if 'RNAup' in self.RNAup and len(inputseq) > 2000:
                self.task_queue.task_done()
                if get_fold:
                    self.result_queue.put(next_task[2:]+[0, '', ''])         
                else:
                    self.result_queue.put(next_task[1:]+[0])
                continue
            if inputseq[-1] != '\n': inputseq += '\n'
            #remove gaps
            inputseq=inputseq.replace('-','').replace('.','')

            try:
                self.p.stdin.write(inputseq)
                self.p.stdin.flush()

            except IOError:
#                self.task_queue.put(next_task)
                self.task_queue.task_done()
                if get_fold:
                    self.result_queue.put(next_task[2:]+[0, '', ''])
                else:
                    self.result_queue.put(next_task[1:]+[0])
                break
            counter+=1
            nrgline = self.p.stdout.readline().rstrip()
#            sys.stderr.write('%s%s'%(inputseq, nrgline))
            bulk=1
            if 'RNAduplex' in self.RNAup:
                bulk=0

            for i in range(bulk):
                fold_line = self.p.stdout.readline()
            if ('inf' in nrgline):
                self.p.stdin.write('@\n')
                self.p = Popen(self.RNAup, shell=True, stdin=PIPE,
                               stdout=PIPE, stderr=PIPE, close_fds=True)
            try:
                nrg = float(
                    nrgline.split(None,1)[1].split('(')[1].split(')')[0].split()[0])
            except IndexError:
                #self.task_queue.put(next_task)
                self.task_queue.task_done()
                if get_fold:
                    self.result_queue.put(next_task[2:]+[0, '', ''])
                else:
                    self.result_queue.put(next_task[1:]+[0])
                #kill this
#                try:
#                    self.p.stdin.write('@\n')
#                except:
#                    pass
#                break
            else:
                self.task_queue.task_done()
                if get_fold:
                    self.result_queue.put(next_task[2:]+[nrg, nrgline, fold_line])
                else:
                    self.result_queue.put(next_task[1:]+[nrg])
            if counter>1000:
                self.p.stdin.write('@\n')
                self.p = Popen(self.RNAup, shell=True, stdin=PIPE,
                               stdout=PIPE, stderr=PIPE, close_fds=True)
                counter=0
#end of class Runner

class RNAupTarPred():
    '''
    implement the TarPred class with RNAup. By default the RNAup 1.8 is used
    with the -b option meaning the energy to open the sRNA is also calculated
    '''
    def __init__(self, cmd='RNAup -Xp -w 25 -b -o',servers=None):
        self.cmd=cmd
        self.servers = servers
        self.initRunners()


    def initRunners(self):
        """
        Initialises the Runners, called once for each object, can be transferred between
        objects
        """
        self.tasks = Queue()
        self.results = Queue()
        if isinstance(self.servers,list):
            self.runners  = [Runner(self.tasks,self.results,
                               'ssh '+server+' nice +10 '+self.cmd)
                        for server in self.servers]
        else:
            try:
#                if isinstance(self.servers,int) or self.servers is None:
                nthreads = int(self.servers)
            except TypeError:
                nthreads = None
            if nthreads is None : nthreads=1
            self.runners = [Runner(self.tasks,self.results,self.cmd) for _ in xrange(nthreads)]
        for w in self.runners:
            w.start()

    def killRunners(self):
        """
        Should be called after the object is used to kill all the threads
        """
        #kill all runners:
        for runner in self.runners:
            self.tasks.put(None)

    def __del__(self ):
        """
        kill the runners
        """
        self.killRunners()



    def scoreall(self,sseq,mseqs):
        '''
        overrides the scoeall of TarPredTop in a way that will run all
        the genome threaded
        '''
        #initialize the runners

        numjobs  = 0
        tarscore = {}
        for gene in mseqs:
            if len(mseqs[gene]) and len(sseq):
                self.tasks.put([str(sseq)+'\n'+str(mseqs[gene])+'\n',gene])
                numjobs+=1
            else:
                tarscore[gene] = 0
        while numjobs:
            result = self.results.get()
            numjobs -= 1
            [gene,nrg] = result
            tarscore[gene] =- nrg
        return tarscore    
        


    def score_pairs(self,seq_pairs, get_fold=False):
        '''
        overrides the scoeall of TarPredTop in a way that will run all
        the genome threaded
        '''
        #initialize the runners

        numjobs  = 0
        tarscore = {}
        folds = {}
        for pname, seqs in seq_pairs.items():
            if len(seqs[0]) and len(seqs[1]):
                if get_fold:
                    self.tasks.put(
                        [str(seqs[0])+'\n'+str(seqs[1])+'\n','get_fold', pname])
                else:
                    self.tasks.put([str(seqs[0])+'\n'+str(seqs[1])+'\n',pname])
                numjobs+=1
            else:
                tarscore[pname] = 0
        while numjobs:
            result = self.results.get()
            numjobs -= 1
            if get_fold:
                [gene, nrg, nrgline, foldline] = result
                folds[gene] = nrgline+'\n'+foldline
            else:
                [gene,nrg] = result
            tarscore[gene] =- nrg
        if get_fold:
            return tarscore, folds
        return tarscore    
        

        
    
