#!/usr/bin/env python
# coding: utf-8

# In[93]:


#نادين اشرف كمال المغربي
#هبه الله محمدممدوح
#ريم وليد علي الشناوي
#مصطفي عادل ذكريا حسن
#مريم ايمن محمد السيد


# In[94]:


import numpy as np
import pandas as pd
import csv
import itertools
from itertools import permutations


# In[95]:


BlosumMatrix = pd.read_csv("blosum62.csv")
BlosumMatrix= BlosumMatrix.iloc[: , :-1]
BlosumMatrix= BlosumMatrix.iloc[: , :-1]
#removethe last two rows
BlosumMatrix = BlosumMatrix.iloc[:-2]
#remove the last two columns
#BlosumMatrix


# In[96]:


newlist=[]
def DeleteLowCompXregion(query,minumim,maximum,idx):
    if maximum + 1 > len(query) or maximum+ 1 == len(query):
        if minumim!=len(query):
            newlist.insert(idx,query[minumim])
            DeleteLowCompXregion(query, minumim + 1, maximum + 1, query + 1)
    elif query[minumim] == query[maximum] and query[minumim + 1] == query[maximum + 1]:
         DeleteLowCompXregion(query,minumim+2,maximum+2,idx)
    elif query[maximum]==query[maximum-2] :
         DeleteLowCompXregion(query, minumim + 2, maximum + 2, idx)
    elif query[minumim]==query[maximum] and query[minumim+1]!=query[maximum+1]:
         DeleteLowCompXregion(query, minumim + 1, maximum + 2, idx)

    else:
        newlist.insert(idx, query[minumim])

        removeRebeated(query, minumim + 1, maximum + 1,idx+1)
    #print(newlist)
    return newlist


# In[97]:


def BlosumDictionary(matrix):#function that takes the BlosumMatrix Dataframe and convert it to dictionary
    AminoAcids_Scores=BlosumMatrix.values.tolist()#list of all the scores.
    charVecs = []
    for ch in range(0,2):#for loop that append the aminoacid letters twice in charVecs list
        charVecs.append(list(matrix.index))
    #split them in 2 seprated lists to take all the possible combinations of 2 aminoacidletters
    AminoAcidLetters1 = charVecs[0]
    AminoAcidLetters2 = charVecs[1]
    AminoAcidCombined= list(itertools.product(AminoAcidLetters1,AminoAcidLetters2))
    #to make the amino acid list of list,each list contains the combination of one letter with all the letters
    #to match each index to it to amino acid scores.
    K = 22
    AminoAcidKeys= []
    for idx in range(0, K):#to make the amino acid list of list,each list contains the possible amino acid letters
        AminoAcidKeys.append( AminoAcidCombined[idx::K])
    BlosumDictionary=dict(zip(AminoAcidKeys[0],AminoAcids_Scores[0]))#initiliaze the
    for i in range(1,22):
        IteratedDictionary= dict(zip(AminoAcidKeys[i],AminoAcids_Scores[i]))#convert 2 lists to dictionary
        BlosumDictionary.update(IteratedDictionary)#append them
    return BlosumDictionary
    


# In[98]:


s=BlosumDictionary(BlosumMatrix)
s[("A","B")]


# In[99]:


def calcScore(word1,word2):#function takes 2 words and calculate their total score
    Blosum_Dict=BlosumDictionary(BlosumMatrix)
    Total_score=0;
    for i in range(0,len(word1)):
        SingleScore=Blosum_Dict[(word1[i],word2[i])]
        Total_score=Total_score+SingleScore
    return Total_score
     


# In[100]:


def SlicingTheWord(slicingNumber,query):#function that takes the word length to output(the slicing words from the query)
    word_list = []
    for wordIndex in range(0,len(query)-slicingNumber+1):#subtract the word length and add 1 to it will give the range of the word
        word_list.append(query[wordIndex:(wordIndex+slicingNumber)])#each word start from the wordIndex value#for example 0
    return word_list                                                #and then 0+slicing number to get a one word


# In[101]:


def getScore(letter1,letter2):#function compute score o 2 aminoacid letters
    Blosum_Dict=BlosumDictionary(BlosumMatrix)
    score=Blosum_Dict[(letter1,letter2)]
    return score
    


# In[102]:


def wordPossibilites(word,matrix): #function that takes a one word and the blosum matrix to get all the possible words of it   
    wordList = []
    for wordIndex in range(0,len(word)):#append each letter of a word into list called wordList
        wordList.append(word[wordIndex])
    #print(wordList)
    charVecs = []    
    for ch in wordList:
        charVecs.append(list(matrix.index))#get all aminoacid letters in list
    #print(charVecs)
    listofPossibility=[]#lest that will have all the possible words of a word
    for i in range(0,len(word)):
        for j in range(0,len(charVecs[0])):
            list1=[]
            list1[:0]=word#take the word into a list
            list1[i]=charVecs[0][j]#iterate through aminoacid sequences and put it instead of an aminoacid letter in the word
            listofPossibility.append(list1)#to have a list contains list of each possible combinantion
    return listofPossibility     


# In[103]:


def GetSeeds(word,wordPossibilites,threshold,wordlength):#function th
    Blosum_Dict=BlosumDictionary(BlosumMatrix)
    seeds=[]
    #seeds=set()
    for j in range(0,len(wordPossibilites)):#loop through all the possible words 
        Total_score=0
        SingleAminoAcid=wordPossibilites[j]#a one possible word
        for i in range(0,wordlength):#when this for loop ends the outer will be incremented to get the j'th word
            SingleScore=Blosum_Dict[(word[i],SingleAminoAcid[i])]#get the score of each letter in the possible word
            Total_score=Total_score+SingleScore      #with the corresponding letter in the original wordword
            if(Total_score>=threshold):#if total score exceeds the threshold
                seeds.append(SingleAminoAcid)#append it in the seeds list
    FinalSeeds = []
    for i in seeds:#if the seed is repeated do not take it as it comes from the same word
        if i not in FinalSeeds:
            FinalSeeds.append(i)
    
    return FinalSeeds


# In[104]:


def listToString(s):#function convert list to string 
    
    # initialize an empty string
    str1 = "" 
    
    # traverse in the string  
    for ele in s: 
        str1 += ele  
    
    # return string  
    return str1 


# In[105]:



def getIndex(string, substring):#function that get the index of substring in string
    index = 0
    if substring in string:
        c = substring[0]
        for ch in string:
            if ch == c:
                if string[index:index+len(substring)] == substring:
                    return index

            index += 1

    return -1


# In[106]:


def ExtendAlignment(query,A_One_Seed, databasesequence,Theword,QueryIndex1,IndexOfTheSeed_InDatabase,score,Hsspthreshold,wordLength):
    QueryIndex2=QueryIndex1+(wordLength-1)
    indexFound2=IndexOfTheSeed_InDatabase+(wordLength-1)
    hsspindex=0
    HSSPSequence=[]
    HSSPScores=[]
    HSSPScores.append(score)#append the score of the word with the query to add to it the right and left scores
    ExtendRightScore=0
    ExtendLeftScore=0
    TotalScore=0
    while True:
        RightExist=False
        LeftExist=False
        if QueryIndex1!=0:#if the seed word is not the first letter of the query
            QueryIndex1decrement=QueryIndex1-1
            queryLetter=query[QueryIndex1decrement] #take the letter before it for extension
            if IndexOfTheSeed_InDatabase!=0:#and also if the seed in the database is not the first letter
                IndexOfTheSeed_InDatabasedecrement=IndexOfTheSeed_InDatabase-1
                databaseLetter=databasesequence[IndexOfTheSeed_InDatabasedecrement]#take the letter before it for extension
                ExtendLeftScore=getScore(queryLetter,databaseLetter)
                LeftExist=True#if both found so we will extend left so left is true and we eill decrement both
                QueryIndex1=QueryIndex1-1
                IndexOfTheSeed_InDatabase=IndexOfTheSeed_InDatabase-1    
        if QueryIndex2!=len(query)-1:#if the last index of the seed word in the query is not the last index in the query
            QueryIndex2increment=QueryIndex2+1
            queryLetter2=query[QueryIndex2increment] #get the letter after it for extension
            if indexFound2!=len(databasesequence)-1:#same for database if the seed is not found at the end of sequence 
                Indexfound2increment=indexFound2+1
                databaseLetter2=databasesequence[Indexfound2increment]#take letter after it
                ExtendRightScore=getScore(queryLetter2,databaseLetter2)
                RightExist=True#extend to right
                QueryIndex2=QueryIndex2+1
                indexFound2=indexFound2+1
        if(LeftExist==True and RightExist==True):#if both are true
            TotalScore3=ExtendLeftScore+ExtendRightScore#get score of both left and right
            Hsspseq=databasesequence[IndexOfTheSeed_InDatabase:indexFound2+1]#get a string of the sequence in db that include the left and right letters
            for i in range(0,len(HSSPScores)):#for loop will add the score of left and right to the previous scores
                Totalscore=TotalScore3+HSSPScores[i]
            HSSPSequence.append(Hsspseq)#list for hssp sequences,append to it the sequence
            HSSPScores.append(Totalscore)#list for hssp scores ,append to it the scores
            Diff=HSSPScores[-1]-HSSPScores[-2]#calculate the difference between the last 2 scores
            Diff=abs(Diff)
            if(Diff>=Hsspthreshold and(HSSPScores[-1]<HSSPScores[-2])):#if difference is greater and the last score is less than the previous
                X=len(HSSPScores)-HSSPScores[::-1].index(max(HSSPScores))-1#take the maximum index of hssp scores from the end
                maxpos2 = X-1#get the maximum position of hssp sequence which lies an index before the score ,because the score is initiliazed by the seed score with the word in the query
                print("Highest score segment pair value",HSSPScores[X]," ")
                if(maxpos2==-1):
                    print("Highest score segment pair sequence ",A_One_Seed)
                else:
                    print("Highest score segment pair sequence ",HSSPSequence[maxpos2])
                break
               #break since the difference exceeds the threshold
   
        elif(RightExist==True):
            TotalScore2=ExtendRightScore
            Hsspseq2=databasesequence[IndexOfTheSeed_InDatabase:indexFound2+1]
            for i in range(0,len(HSSPScores)):
                Totalscore=TotalScore2+HSSPScores[i]
            HSSPSequence.append(Hsspseq2)
            HSSPScores.append(Totalscore)
            Diff=HSSPScores[-1]-HSSPScores[-2]
            Diff=abs(Diff) 
            if(Diff>=Hsspthreshold and(HSSPScores[-1]<HSSPScores[-2])):
                maxvalue=max(HSSPScores)
                maxpos=HSSPScores.index(maxvalue)
                X=len(HSSPScores)-HSSPScores[::-1].index(max(HSSPScores))-1
                maxpos2 = X-1
                print("Highest score segment pair value",HSSPScores[X]," ")
                if(maxpos2==-1):
                    print("Highest score segment pair sequence ",A_One_Seed)
                else:
                    print("Highest score segment pair sequence ",HSSPSequence[maxpos2])
                break
        elif(LeftExist==True):
            TotalScore1=ExtendLeftScore
            Hsspseq3=databasesequence[IndexOfTheSeed_InDatabase:indexFound2+1]
            for i in range(0,len(HSSPScores)):
                Totalscore=TotalScore1+HSSPScores[i]
            HSSPSequence.append(Hsspseq3)
            HSSPScores.append(Totalscore)
            Diff=HSSPScores[-1]-HSSPScores[-2]
            Diff=abs(Diff)
            if(Diff>=Hsspthreshold and(HSSPScores[-1]<HSSPScores[-2])):
                maxvalue=max(HSSPScores)
                maxpos=HSSPScores.index(maxvalue)
                X=len(HSSPScores)-HSSPScores[::-1].index(max(HSSPScores))-1
                maxpos2 = X-1
                print("Highest score segment pair value",HSSPScores[X]," ")
                if(maxpos2==-1):
                    print("Highest score segment pair sequence ",A_One_Seed)
                else:
                    print("Highest score segment pair sequence ",HSSPSequence[maxpos2])
                break
        else:#if right and left is false#the sequences /query ends
            maxvalue=max(HSSPScores)
            maxpos=HSSPScores.index(maxvalue)
            X=len(HSSPScores)-HSSPScores[::-1].index(max(HSSPScores))-1
            maxpos2 = X-1
            print("Highest score segment pair value",HSSPScores[X]," ")
            if(maxpos2==-1):
                print("Highest score segment pair sequence ",A_One_Seed)
            else:
                print("Highest score segment pair sequence ",HSSPSequence[maxpos2])
            break
        


# In[107]:


def BlastAlgorithm(BlosumMatrix):
    df=pd.read_table("C:\\Users\\kc\\Desktop\\seed.txt")
    Dbseq=df.values.tolist()#convert the textfile to list of list each list is a DB sequence
    #Dbseq=[["ALNATYGTYAMNBAZ"],["AMNGTYGAYBCD"],["AMNBAZGNYBAZ"]]
    Query = input("Enter your query : ")
    Query=Query.upper()
    query=""
    char=['A','C','D','B','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
    for i in range (0,len(Query)):
        if(Query[i] not in char):
            Query = input("Enter valid query : ")
            Query=Query.upper()
    #print(Query)
    wordLengthInput=input("Enter the word length: ")
    wordLength=int(wordLengthInput);
    threshold2=input("Enter the threshold")
    threshold=int(threshold2)
    hsspthreshold2 = input("Enter your hsspth : ")
    hsspthreshold=int(hsspthreshold2)
    #newsequence=[]
    newsequence=DeleteLowCompXregion(Query,0,2,0)
    #print("outPut of step 1")
    #print(newsequence)
    query=listToString(newsequence)
    #print("SEQ")
    #print( query)
    word=SlicingTheWord(wordLength,query)
    Finalseeds=[]
    queryexists=False
    for i in range(0,len(word)):
        SlicedWords=wordPossibilites(word[i],BlosumMatrix)
        listofseeds=GetSeeds(word[i],SlicedWords,threshold,wordLength)#this function return the seeds that exceed the thresholf for  A ONE WORD ONLY.
        Finalseeds.append(listofseeds)#so i am appending all of the seeds for all the words in another list
    for queryindex in range (0,len(Finalseeds)):#iterate through all of the seeds(list containg list of seeds for each word)
        for k in range(0,len(Finalseeds[queryindex])):#iterating though a seeds for word number j
            score=0
            A_One_Seed=Finalseeds[queryindex][k]
            A_One_Seed=listToString(A_One_Seed)#get a one seed in a string
            for Databaseseq in range(0,len(Dbseq)):#loop through the whole database sequences for 1 seed
                databasesequence=listToString(Dbseq[Databaseseq])#convert DB sequence to list
                IndexOfTheSeed_InDatabase=getIndex(databasesequence,A_One_Seed)#get the index of first letter of the seed in the Database
                if(IndexOfTheSeed_InDatabase!=-1):#if the seed exists
                    queryexists=True
                    print("The seed ",A_One_Seed," is found in database sequence number: ", Databaseseq+1)
                    score=calcScore(A_One_Seed,word[queryindex])
                    ExtendAlignment(query,A_One_Seed,databasesequence,word[queryindex],queryindex,IndexOfTheSeed_InDatabase,score,hsspthreshold,wordLength)
                   #call extend alignment method.
    if(queryexists==False):
        print("Query not found in the database")


# In[108]:


BlastAlgorithm(BlosumMatrix)


# # 
