
# coding: utf-8

# In[1]:


from sklearn.ensemble import RandomForestClassifier
import sys, re
import numpy as np
from sklearn import metrics
from sklearn.model_selection import train_test_split
import joblib
import sys

overlap_model = sys.argv[4]
nonoverlap_model = sys.argv[5]
clf_overlap = joblib.load(overlap_model)
clf_nonoverlap = joblib.load(nonoverlap_model)



def split_query_file_to_features(query_file):
        X = []
        with open(query_file, 'r') as f:
                for line in f:
                        sp = line.strip().split('\t')
                        if sp[0] == 'id': continue
                        #X.append(np.array(list(map(float, sp[1:len(sp)]))))
                        X.append(np.array([float(s) for s in sp[1:len(sp)]]))
        return np.array(X)


path = sys.argv[1]
save_path = path
overlap_pass = sys.argv[2]
query_file = path + overlap_pass
X_query = split_query_file_to_features(query_file)
if X_query.shape[0] == 0:
	open(query_file + ".RFpred.csv", 'w').close()
else:
	y_query_predict=clf_overlap.predict_proba(X_query)[:,1]
	save_filename = query_file + ".RFpred.csv"
	np.savetxt(save_filename, y_query_predict, delimiter=",")




def split_query_file_to_features(query_file):
        X = []
        with open(query_file, 'r') as f:
                for line in f:
                        sp = line.strip().split('\t')
                        if sp[0] == 'id': continue
                        #X.append( np.array(map(float, sp[1:len(sp)])) )
                        X.append(np.array([float(s) for s in sp[1:len(sp)]]))
        return np.array(X)




nonoverlap_pass = sys.argv[3]
query_file = path + nonoverlap_pass
X_query = split_query_file_to_features(query_file)

if X_query.shape[0] == 0:
	open(query_file + ".RFpred.csv", 'w').close()
else:
	y_query_predict=clf_nonoverlap.predict_proba(X_query)[:,1]
	save_filename = query_file + ".RFpred.csv"
	np.savetxt(save_filename, y_query_predict, delimiter=",")





