import numpy as np
import warnings
from math import log, exp
from pprint import PrettyPrinter


#np.random.choice
def main():
    states = [
        State("A", {"A": 0.98,
                    "H": 0.01,
                    "B": 0.01}, {"a": 0.9,
                                 "b": 0.1}),
        State("H", {"A": 0.01,
                    "H": 0.98,
                    "B": 0.01}, {"a": 0.5,
                                 "b": 0.5}),
        State("B", {"A": 0.01,
                    "H": 0.01,
                    "B": 0.98}, {"a": 0.1,
                                 "b": 0.9})
    ]

    hmm1 = HMM(states, {"A": 1 / 3.0, "B": 1 / 3.0, "H": 1 / 3.0})
    obs, sts = hmm1.simulate(3000)

    hmm2 = HMM()
    hmm2.trains(obs, sts)
    for key, state in hmm2.states.items():
        print(state.name)
        print({key: exp(val) for key, val in state.emission.items()})
        print({key: exp(val) for key, val in state.transition.items()})
    
    snps = list(zip(range(len(obs)),obs))
    # print("".join(obs)+"\n"+"".join(sts)+"\n"+"".join(hmm2.viterbi(snps)))
    print(snps)
    print(hmm2.viterbi(snps))

printer = PrettyPrinter(indent=4, width=80)

class State(object):
    """states for use in _HMM"""

    def __init__(self,
                 name,
                 transition={},
                 emission={},
                 observations=float('inf')):
        self.name = name
        self.transition = { 
            key: log(value) # {state_name:P(t)}
            for key, value in transition.items()
        } 
        self.emission = { 
            key: log(value) # {observation_name:P(e)}
            for key, value in emission.items()
        }
        self.observations = observations if (self.emission
                                             or self.transition) else 0
        self.in_state = False

    def sample_obs(self):
        obs, probs = zip(*self.emission.items())
        probs = [exp(p) for p in probs]
        return np.random.choice(obs, p=probs)

    def sample_state(self):
        states, probs = zip(*self.transition.items())
        probs = [exp(p) for p in probs]
        return np.random.choice(states, p=probs)
        
    def p(self, obs):
        return self.emission[obs] if obs in self.emission else 0
    def t(self, state):
        return self.transition[state] if state in self.transition else 0

    def train(self, obs=False, trans=False):
        if (obs and trans): raise Exception("Train observation OR transition")
        if self.observations != float('inf'):
            # update transition probability if not first observation in a sequence
            # this means a transition can only be trained (setting in_state to false)
            # after at least one observation
            if self.in_state:
                if trans: self.in_state = False
                #calculate iterative mean adjustments for current count
                pre = (self.observations - 1) / self.observations
                mean_multiply = log(pre) if (pre != 0) else 0
                mean_add = log(1 / self.observations)

                # update transition probabilities
                for key in self.transition:
                    self.transition[key] += mean_multiply
                trans_name = trans if trans else self.name
                if trans_name not in self.transition:
                    self.transition[trans_name] = mean_add
                else:
                    self.transition[trans_name] = \
                        np.logaddexp(mean_add, self.transition[trans_name])
            # then, for emmisions
            if obs:
                # note we are now in the state.
                self.in_state = True
                # iterate observation count
                self.observations += 1
                #calculate iterative mean adjustments for count
                pre = (self.observations - 1) / self.observations
                mean_multiply = log(pre) if (pre != 0) else 0
                mean_add = log(1 / self.observations)

                #emission observation
                for key in self.emission:
                    self.emission[key] += mean_multiply
                if obs not in self.emission:
                    self.emission[obs] = mean_add
                    # print(obs,":",exp(self.emission[obs]))
                else:
                    self.emission[obs] = \
                        np.logaddexp(mean_add, self.emission[obs])
                    # print(obs,":",exp(self.emission[obs]))
        else:
            warnings.warn("Cannot train without a known observation count!")


class HMM(object):
    """A basic hidden markov model class. Uses log(prob) to reduce floating point precision problems."""

    def __init__(self, states=[], priors={}):
        self.priors = {key: log(value) for key, value in priors.items()}
        self.states = {state.name: state for state in states}
        self.observable = set(obs for state in self.states.values() for obs in state.emission)
        self.last_trained_state = None
        
    def dump(self):
        return {"states": {
                    state.name: {
                        "transition":{
                            key:exp(val) 
                            for key,val in state.transition.items()
                        },
                        "emission":{
                            key:exp(val) 
                            for key,val in state.emission.items()
                        },
                    } 
                    for state in self.states.values()
                },
                "priors":{
                    state.name:exp(self.prior(state.name))
                    for state in self.states.values()
                }}
    
    def load(self,config_dict):
        states = [State(name,state["transition"],state["emission"]) 
                  for name,state in config_dict["states"].items()]
        self.__init__(states=states,priors=config_dict["priors"])
    
    def prior(self,state):
        if self.priors:
            return self.priors[state] if state in self.priors else 0
        else:
            return log(1/len(self.states))
    
    def simulate(self, length):
        states, probs = zip(*self.priors.items())
        probs = [exp(p) for p in probs]
        initial_state = self.states[np.random.choice(states, p=probs)]
        initial_obs = initial_state.sample_obs()
        obs = [initial_obs]
        sts = [initial_state]
        for i in range(length - 1):
            next_state = self.states[sts[-1].sample_state()]
            next_obs = next_state.sample_obs()
            obs.append(next_obs)
            sts.append(next_state)
        sts = [st.name for st in sts]
        return obs, sts
    
    def train_begin(self):
        self.last_trained_state = None
    
    def train(self, ob, st):
        if st not in self.states:
            self.states[st] = State(st)
        curr_state = self.states[st]
        if self.last_trained_state and self.last_trained_state != curr_state:
            self.last_trained_state.train(trans=st)
        curr_state.train(ob)
        self.observable.add(ob)
        self.last_trained_state = curr_state
    
    def train_end(self):
        self.last_trained_state = None
        
    def trains(self, obs, sts):
        train_begin()
        for ob, st in zip(obs, sts):
            train(obs,sts)
        train_end()

    def get_most_likely_state(self,observations):
        '''Calculates the probability a list of SNPs is of each state then returns the most probable state.'''
        return max((sum(state.p(obs) for obs in observations), state.name) for state in self.states.values())[1]
    
    def viterbi(self,observations):
        """"""
        vtrb = []
        trace_graph = {}
    
        snpList = [obs for obs in observations if obs in self.observable]
        
        enum_states = list(enumerate(self.states.values()))
    
        for t,snp in enumerate(snpList):
            obs = snp
            vtrb.append([])
            if t==0:
                for i,state in enum_states:
                    vtrb[0].append(self.prior(state.name) + state.p(obs))
                    trace_graph[(0,i)] = None
            else:
                for i,state in enum_states:
                    max_prob,max_prev = max((vtrb[-2][p] + prev_state.t(state.name) + state.p(obs), p) for p,prev_state in enum_states)
                    vtrb[t].append(max_prob)
                    trace_graph[(t,i)] = (t-1,max_prev)
    
        max_final,last = max((vtrb[-1][i],(len(vtrb)-1,i)) for i,state in enum_states)
        state_dets = []
        prev = last

        while prev!=None:
            state_dets.append((prev[0],enum_states[prev[1]][1].name))
            prev = trace_graph[prev]
        state_dets[:] = state_dets[::-1]

        return state_dets
    
    def __str__(self):
        return printer.pformat(self.dump())
    
    def copy():
        new_hmm = HMM()
        new_hmm.load(self.dump())
        return new_hmm


if __name__ == '__main__':
    main()
