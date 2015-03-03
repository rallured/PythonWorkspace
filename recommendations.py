from numpy import *
import scipy.stats as stat
import pdb

# A dictionary of movie critics and their ratings of a small
# set of movies
critics={'Lisa Rose': {'Lady in the Water': 2.5, 'Snakes on a Plane': 3.5,
 'Just My Luck': 3.0, 'Superman Returns': 3.5, 'You, Me and Dupree': 2.5,
 'The Night Listener': 3.0},
'Gene Seymour': {'Lady in the Water': 3.0, 'Snakes on a Plane': 3.5,
 'Just My Luck': 1.5, 'Superman Returns': 5.0, 'The Night Listener': 3.0,
 'You, Me and Dupree': 3.5},
'Michael Phillips': {'Lady in the Water': 2.5, 'Snakes on a Plane': 3.0,
 'Superman Returns': 3.5, 'The Night Listener': 4.0},
'Claudia Puig': {'Snakes on a Plane': 3.5, 'Just My Luck': 3.0,
 'The Night Listener': 4.5, 'Superman Returns': 4.0,
 'You, Me and Dupree': 2.5},
'Mick LaSalle': {'Lady in the Water': 3.0, 'Snakes on a Plane': 4.0,
 'Just My Luck': 2.0, 'Superman Returns': 3.0, 'The Night Listener': 3.0,
 'You, Me and Dupree': 2.0},
'Jack Matthews': {'Lady in the Water': 3.0, 'Snakes on a Plane': 4.0,
 'The Night Listener': 3.0, 'Superman Returns': 5.0, 'You, Me and Dupree': 3.5},
'Toby': {'Snakes on a Plane':4.5,'You, Me and Dupree':1.0,'Superman Returns':4.0}}

#Get Euclidean distance metric
def euclidean(prefs,person1,person2):
    dist = 0.
    for item in prefs[person1]:
        if item in prefs[person2]:
            dist += sqrt(prefs[person1][item]**2 + prefs[person2][item]**2)

    return 1/(1+dist)

#Get Pearson correlation coefficient
def pearson(prefs,person1,person2):
    x = []
    y = []
    for item in prefs[person1]:
        if item in prefs[person2]:
            x.append(prefs[person1][item])
            y.append(prefs[person2][item])
    return stat.pearsonr(x,y)[0]

#Get ordered list of closest matches
def matches(prefs,person,n=5,similarity=pearson):
    scores = [(similarity(prefs,person,other),other) for other in prefs\
              if other != person]
    scores.sort()
    scores.reverse()
    
    return scores

#Return list of movies along with their ranking score
#Ranking score is weighted average by similarity score
def getRec(prefs,person):
    totals = {} #Store total score
    simSums = {} #Store similarity sum (weight)
    for other in prefs:
        #Don't compare to myself
        if other is not person:
            sim = pearson(prefs,person,other)
            if logical_or(sim < 0.,isnan(sim)):
                continue
            for item in prefs[other]:
                if item not in prefs[person]:
                    totals.setdefault(item,0)
                    totals[item]+=prefs[other][item]*sim
                    simSums.setdefault(item,0)
                    simSums[item]+=sim
    recs = [(totals[item]/simSums[item],item) for item in totals]
    recs.sort()
    recs.reverse()
    return recs

#Swap people with movies
def transform(prefs):
    result = {}
    for person in prefs:
        for item in prefs[person]:
            result.setdefault(item,{})
            result[item][person] = prefs[person][item]
    return result

def calcSimilarItems(prefs,n=10):
    result = {}

    itemPrefs = transform(prefs)
    c=0
    for item in itemPrefs:
        c+=1
        if c%100==0: print "%d / %d" % (c,len(itemPrefs))
        scores = matches(itemPrefs,item,n=n,similarity=pearson)
        result[item] = scores
    return result

def getRecItems(prefs,itemMatch,user):
    userRatings = prefs[user]
    scores = {}
    totalSim = {}

    for (item,rating) in userRatings.items():
        for (similarity,item2) in itemMatch[item]:
            if item2 in userRatings: continue
            scores.setdefault(item2,0)
            scores[item2]+=similarity*rating
            totalSim.setdefault(item2,0)
            totalSim[item2]+=similarity
            
    rankings = [(score/totalSim[item],item) for item,score in scores.items()]
    rankings.sort()
    rankings.reverse()
    return rankings

def loadMovieLens(path='/Users/ryanallured/Collective/ml-100k/'):
    movies = {}
    for line in open(path+'/u.item'):
        (id,title)=line.split('|')[0:2]
        movies[id] = title

    prefs = {}
    for line in open(path+'/u.data'):
        (user,movieid,rating,ts)=line.split('\t')
        prefs.setdefault(user,{})
        prefs[user][movies[movieid]]=float(rating)
    return prefs
