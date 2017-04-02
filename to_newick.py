#!/usr/bin/python2.7

from ete2 import Tree, TreeStyle


seq = '((d:123,p:138):-38,((m:80,k:84):23,(((h:192,u:265):29,(v:52,c:58):15):-2,((l:44,(r:328,e:510):39):3,(b:37,j:37):-3):2):34):-1,((((a:196,i:169):28,((o:607,q:774):15,(n:104,w:101):28):23):8,(g:87,f:90):32):-20,(s:73,t:68):-15):-20);'
t = Tree(seq)
ts = TreeStyle()
ts.show_leaf_name = True
ts.branch_vertical_margin = 10
ts.show_scale = False
ts.show_branch_length = True
t.show(tree_style=ts)