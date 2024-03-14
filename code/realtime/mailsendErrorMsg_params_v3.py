import re
import sys
from email.mime.base import MIMEBase
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
import smtplib, ssl
from email import encoders

import subprocess

import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from io import StringIO
from IPython.display import SVG
import pydot
from treelib import Node, Tree
from fpdf import FPDF

import unicodedata

def unicode_normalize(s):
       return unicodedata.normalize('NFKD', s).encode('ascii', 'ignore')
# python front_end/sendMail_v2.py niladri.das.2010@gmail.com /data2/nmr/our_algo_final/protein/2m4k_param_2m4k/sendfiles/sendfiles_list.txt 2m4k

def successOrfail_one(filename, srch_pattern, success_msg, fail_msg, abscent_msg):
       '''

       '''
       regex = re.compile(srch_pattern)       

       ret_list = []

       with open(filename, 'r') as fp:            
            #ret_list.append(srch_pattern)
            for line in fp:
                  k = regex.findall(line)
                  position = False
                  if str(k) != '[]':
                        position = True
                        if "success" in line.lower():                                                              
                              ret_list.append(success_msg)
                        elif "error" in line.lower():
                              ret_list.append(fail_msg)
                        else:
                              ret_list.append(abscent_msg)
       if not ret_list:
            ret_list = [abscent_msg] 

       return ret_list

def successOrfail_multi(filename, srch_pattern, success_msg, fail_msg, abscent_msg):
       '''

       '''
       regex = re.compile(srch_pattern)

       ret_list = []

       c=0
       with open(filename, 'r') as fp:
            #ret_list.append(srch_pattern)
            for line in fp:
                  k = regex.findall(line)
                  position = False
                  if str(k) != '[]':
                        position = True
                        c = c+1
                        if "success" in line.lower():
                               ret_list.append('%s for %d'%(success_msg,c))
                        elif "error" in line.lower():
                               ret_list.append('%s for %d'%(fail_msg,c))                               
                        else:
                               ret_list.append('%s for %d'%(abscent_msg,c))                               
       if not ret_list:
            ret_list = [abscent_msg]

       return ret_list

def returnMatchingFIlenames(filename, srch_pattern):
      '''

      '''
      ret_list = []

      regex = re.compile(srch_pattern)
      with open(filename, 'r') as fp:
            for line in fp:
                  k = regex.findall(line)
                  position = False
                  if str(k) != '[]':
                        position = True
                        parts = line.split()

                        ret_list.append(parts[-1])
      
      return ret_list

def parseMessageFiles_archived(filename, oppfile):
      '''

      '''
      
      parsing_file_pat    = "in parsing the input parameter files"
      parsing_div_frags   = "and individual modeling section"
      parsing_consolidate = "consolidating the divided parts"
      parsing_gaps_model  = "modelled the gaps post processing steps"
      parsing_gaps_post   = "modelled the gaps post processing steps"

      dot_graph = pydot.Dot(graph_type='digraph')
      tree = Tree()
      tree.create_node("DREAM","D0")

      #fp = open(oppfile, 'w')
      #fp.write('\n======================== Parsing parameter file ========================')
      #fp.close()
      pred_stage = pydot.Node('DREAM')
      pred_stage.set_shape('box3d')
      dot_graph.add_node(pred_stage)

      ret_list = successOrfail_one(filename, parsing_file_pat, "Successfully parsed parameter file.", "Error in parsing paramater file.", "Abscent status in parsing parameter.")      
      pred_stage_i = pydot.Node(ret_list[0], margin=0, width=0, height=0)
      if  "Successfully" in ret_list[0]:
            pred_stage_i = pydot.Node('A', style='filled', color='green', margin=0, width=0, height=0)
            #tree.create_node("A: %s"%(ret_list[0]),"A", parent="D0")            
      elif "Error" in ret_list[0]:
            pred_stage_i = pydot.Node('A', style='filled', color='red', margin=0, width=0, height=0)
            #tree.create_node("A: %s"%(ret_list[0]),"A", parent="D0")                        
      elif "Abscent" in ret_list[0]:
            pred_stage_i = pydot.Node('A', style='filled', color='orange', margin=0, width=0, height=0)
            #tree.create_node("A: %s"%(ret_list[0]),"A", parent="D0")
      tree.create_node("A: Parsing parameter file","A", parent="D0")
      pred_stage_i.set_shape('box')
      dot_graph.add_node(pred_stage_i)
      pred_stage_edge = pydot.Edge(pred_stage, pred_stage_i)
      pred_stage_edge.set_label('')
      dot_graph.add_edge(pred_stage_edge)

      #fp = open(oppfile, 'a')
      #fp.write('\n======================== Dividing into fragments and modelling ========================')
      #fp.close()     
      ret_list = successOrfail_one(filename, parsing_div_frags, "Successfully divided into fragments and modelled them.", "Error while dividing into fragments and their modelling.", "Abscent status while dividing into fragments and their modelling.")
      #div_n_calc = pydot.Node(ret_list[0], margin=0, width=0, height=0)
      if  "Successfully" in ret_list[0]:
            div_n_calc = pydot.Node('B', style='filled', color='green', margin=0, width=0, height=0)
            #tree.create_node("B: %s"%(ret_list[0]),"B", parent="D0")
      elif "Error" in ret_list[0]:
            div_n_calc = pydot.Node('B', style='filled', color='red', margin=0, width=0, height=0)
            #tree.create_node("B: %s"%(ret_list[0]),"B", parent="D0")            
      elif "Abscent" in ret_list[0]:
            div_n_calc = pydot.Node('B', style='filled', color='orange', margin=0, width=0, height=0)
            #tree.create_node("B: %s"%(ret_list[0]),"B", parent="D0")      
      tree.create_node("B: dividing into fragments and modelling","B", parent="D0")            
      div_n_calc.set_shape('box')
      dot_graph.add_node(div_n_calc)
      div_n_calc_edge = pydot.Edge(pred_stage_i, div_n_calc)
      div_n_calc_edge.set_label('')
      dot_graph.add_edge(div_n_calc_edge)      

      #fp = open(oppfile, 'a')
      #fp.write('\n======================== Consolidating structures for fragments ========================')
      #fp.close()
      ret_list = successOrfail_one(filename, parsing_consolidate, "Successfully consolidated structure for fragments.", "Error consolited structure for fragments.", "Abscent status while consolidating structure for fragments.")
      #register = pydot.Node(ret_list[0], margin=0, width=0, height=0)
      if  "Successfully" in ret_list[0]:
            register = pydot.Node('C', style='filled', color='green', margin=0, width=0, height=0)
            #tree.create_node("C: %s"%(ret_list[0]),"C", parent="D0")
      elif "Error" in ret_list[0]:
            register = pydot.Node('C', style='filled', color='red', margin=0, width=0, height=0)
            #tree.create_node("C: %s"%(ret_list[0]),"C", parent="D0")            
      elif "Abscent" in ret_list[0]:
            register = pydot.Node('C', style='filled', color='orange', margin=0, width=0, height=0)
            #tree.create_node("C: %s"%(ret_list[0]),"C", parent="D0")         
      tree.create_node("C: consolidating structures for fragments","C", parent="D0")
      register.set_shape('box')
      dot_graph.add_node(register)
      register_edge = pydot.Edge(div_n_calc, register)
      register_edge.set_label('')
      dot_graph.add_edge(register_edge)       

      #fp = open(oppfile, 'a')
      #fp.write('\n======================== Modelling the gaps ========================')
      #fp.close()
      ret_list = successOrfail_one(filename, parsing_gaps_model, "Successfully modelled gaps", "Error while modelling gaps", "Abscent status Error while modelling gaps")
      #gaps = pydot.Node(ret_list[0], margin=0, width=0, height=0)
      if  "Successfully" in ret_list[0]:
            gaps = pydot.Node('D', style='filled', color='green', margin=0, width=0, height=0)
            #tree.create_node("D: %s"%(ret_list[0]),"D", parent="D0")
      elif "Error" in ret_list[0]:
            gaps = pydot.Node('D', style='filled', color='red', margin=0, width=0, height=0)
            #tree.create_node("D: %s"%(ret_list[0]),"D", parent="D0")            
      elif "Abscent" in ret_list[0]:
            gaps = pydot.Node('D', style='filled', color='orange', margin=0, width=0, height=0)
            #tree.create_node("D: %s"%(ret_list[0]),"D", parent="D0")
      tree.create_node("D: Modelling gaps","D", parent="D0")            
      gaps.set_shape('box')
      dot_graph.add_node(gaps)
      gaps_edge = pydot.Edge(register, gaps)
      gaps_edge.set_label('')
      dot_graph.add_edge(gaps_edge)            

      '''
      #fp = open(oppfile, 'a')
      #fp.write('\n======================== Prediction stage ========================')
      #fp.close()
      ret_list = successOrfail_one(filename, oppfile, parsing_gaps_post, "Successfully modelled tha gaps.", "Error while modelling gaps.", "Abscent status Error while modelling gaps.")      
      '''

      parsing_gaps_1_model_i = "generated starting structure generation for gaps"
      parsing_gaps_2_model_i = "in choosing starting structures for gap regions"
      parsing_gaps_3_model_i = "together choosing starting structures for gap regions"
      parsing_gaps_4_model_i = "together chosen structrue for the gaps"
      parsing_gaps_5_model_i = "post processing EM step"
      parsing_gaps_6_model_i = "for water refinement step"
      parsing_gaps_7_model_i = "water refinement step in post processing"

      #fp = open(oppfile, 'a')
      #fp.write('\n======================== Gap correction starting sturcture ========================')
      #fp.close()
      dummy_2 = pydot.Node('1', margin=0, width=0, height=0)
      tree.create_node("E: Starting structure for gaps", "1", parent="D0")
      dummy_2.set_shape('circle')
      dot_graph.add_node(dummy_2)      
      ret_list = successOrfail_multi(filename, parsing_gaps_1_model_i, "Successfully modelled starting structure for gaps", "Error in modelling starting structure for gaps", "Unknown status for modelling starting structure for gaps")
      count = 0
      for ele in ret_list:
            gap_start = pydot.Node('gap correction starting sturcture\n%s'%(ele), margin=0, width=0, height=0)
            count = count+1
            if  "Successfully" in ele:
                  gap_start = pydot.Node("E.1.%d"%(count), style='filled', color='green', margin=0, width=0, height=0)
                  #tree.create_node("E.1.%d: %s"%(count,ele), "E.1.%d"%(count), parent="1")
            elif "Error" in ele:
                  gap_start = pydot.Node("E.1.%d"%(count), style='filled', color='red', margin=0, width=0, height=0)
                  #tree.create_node("E.1.%d: %s"%(count,ele), "E.1.%d"%(count), parent="1")            
            elif "Abscent" in ele:
                  gap_start = pydot.Node("E.1.%d"%(count), style='filled', color='orange', margin=0, width=0, height=0)
                  #tree.create_node("E.1.%d: %s"%(count,ele), "E.1.%d"%(count), parent="1")
            tree.create_node("E.1.%d: Modelling starting structrues for gaps"%(count), "E.1.%d"%(count), parent="1")
            gap_start.set_shape('box')
            dot_graph.add_node(gap_start)
            gap_start_edge = pydot.Edge(gaps, gap_start)
            dot_graph.add_edge(gap_start_edge)
            dummy_2_anchor = pydot.Edge(gap_start, dummy_2)
            dot_graph.add_edge(dummy_2_anchor)
      
      #fp = open(oppfile, 'a')
      #fp.write('\n======================== Generating structures for the gap correct ========================')
      #fp.close()
      dummy_3 = pydot.Node('2', margin=0, width=0, height=0)
      tree.create_node("F: Starting structure for gap correct", "2", parent="D0")
      dummy_3.set_shape('circle')
      dot_graph.add_node(dummy_3)      
      ret_list = successOrfail_multi(filename, parsing_gaps_2_model_i, "Successfully generated starting structures for gaps correct", "Error in generating starting structures for gaps correct", "Unknown status in generating starting structures for gaps correct")
      count = 0
      for ele in ret_list:
            #gap_start2 = pydot.Node('choosen starting structures for gap regions\n%s'%(ele), margin=0, width=0, height=0)
            count = count+1
            if  "Successfully" in ele:
                  gap_start2 = pydot.Node("F.2.%d"%(count), style='filled', color='green', margin=0, width=0, height=0)
                  #tree.create_node("F.2.%d: %s"%(count,ele), "F.2.%d"%(count), parent="2")
            elif "Error" in ele:
                  gap_start2 = pydot.Node("F.2.%d"%(count), style='filled', color='red', margin=0, width=0, height=0)
                  #tree.create_node("F.2.%d: %s"%(count,ele), "F.2.%d"%(count), parent="2")            
            elif "Abscent" in ele:
                  gap_start2 = pydot.Node("F.2.%d"%(count), style='filled', color='orange', margin=0, width=0, height=0)
                  #tree.create_node("F.2.%d: %s"%(count,ele), "F.2.%d"%(count), parent="2")
            tree.create_node("F.2.%d: Structures for gap correct"%(count), "F.2.%d"%(count), parent="2")                              
            gap_start2.set_shape('box')
            dot_graph.add_node(gap_start2)
            gap_start2_edge = pydot.Edge(dummy_2, gap_start2)
            dot_graph.add_edge(gap_start2_edge)
            dummy_3_anchor = pydot.Edge(gap_start2, dummy_3)
            dot_graph.add_edge(dummy_3_anchor)                 

      #fp = open(oppfile, 'a')
      #fp.write('\n======================== Putting together the gaps corrected initial ========================')
      #fp.close()
      dummy_4 = pydot.Node('3', margin=0, width=0, height=0)
      tree.create_node("G: Putting together the chosen structures for gap correct initial stage", "3", parent="D0")
      dummy_4.set_shape('circle')
      dot_graph.add_node(dummy_4)
      ret_list = successOrfail_multi(filename, parsing_gaps_3_model_i, "Successfully put together starting structures for gap correct initial stage", "Error in putting together starting structures for gap regions initial stage", "Unknown status for putting together starting structures for gap regions initial stage")
      count = 0
      for ele in ret_list:
            #put_together_gap = pydot.Node('choosen starting structures for gap regions\n%s'%(ele), margin=0, width=0, height=0)
            count = count+1
            if  "Successfully" in ele:
                  put_together_gap = pydot.Node("G.3.%d"%(count), style='filled', color='green', margin=0, width=0, height=0)
                  #tree.create_node("G.3.%d: %s"%(count,ele), "G.3.%d"%(count), parent="3")
            elif "Error" in ele:
                  put_together_gap = pydot.Node("G.3.%d"%(count), style='filled', color='red', margin=0, width=0, height=0)
                  #tree.create_node("G.3.%d: %s"%(count,ele), "G.3.%d"%(count), parent="3")            
            elif "Abscent" in ele:
                  put_together_gap = pydot.Node("G.3.%d"%(count), style='filled', color='orange', margin=0, width=0, height=0)
                  #tree.create_node("G.3.%d: %s"%(count,ele), "G.3.%d"%(count), parent="3")
            tree.create_node("G.3.%d: Put together starting structure for gap correct initial"%(count), "G.3.%d"%(count), parent="3")
            put_together_gap.set_shape('box')
            dot_graph.add_node(put_together_gap)
            put_together_gap_edge = pydot.Edge(dummy_3, put_together_gap)
            dot_graph.add_edge(put_together_gap_edge)
            dummy_4_anchor = pydot.Edge(put_together_gap, dummy_4)
            dot_graph.add_edge(dummy_4_anchor)  
      
      #fp = open(oppfile, 'a')
      #fp.write('\n======================== Putting together the gaps corrected final ========================')
      #fp.close()
      dummy_5 = pydot.Node('4', margin=0, width=0, height=0)
      tree.create_node("H: Putting together the chosen structures for gap correct final stage", "4", parent="D0")
      dummy_5.set_shape('circle')
      dot_graph.add_node(dummy_5)      
      ret_list = successOrfail_multi(filename, parsing_gaps_4_model_i, "Successfully put together choosen structrue for the gaps final stage", "Error in putting together choosen structrue for the gaps final stage", "Unknown status for together chosen structrue for the gaps final stage")
      count = 0
      for ele in ret_list:
            #put_together_gap2 = pydot.Node('choosen starting structures for gap regions\n%s'%(ele), margin=0, width=0, height=0)
            count = count+1
            if  "Successfully" in ele:
                  put_together_gap2 = pydot.Node("H.4.%d"%(count), style='filled', color='green', margin=0, width=0, height=0)
                  #tree.create_node("H.4.%d: %s"%(count,ele), "H.4.%d"%(count), parent="4")
            elif "Error" in ele:
                  put_together_gap2 = pydot.Node("H.4.%d"%(count), style='filled', color='red', margin=0, width=0, height=0)
                  #tree.create_node("H.4.%d: %s"%(count,ele), "H.4.%d"%(count), parent="4")            
            elif "Abscent" in ele:
                  put_together_gap2 = pydot.Node('H.4', style='filled', color='orange', margin=0, width=0, height=0)
                  #tree.create_node("H.4.%d: %s"%(count,ele), "H.4.%d"%(count), parent="4")
            tree.create_node("H.4.%d: Put together starting structure for gap correct final"%(count), "H.4.%d"%(count), parent="4")                                 
            put_together_gap2.set_shape('box')
            dot_graph.add_node(put_together_gap2)
            put_together_gap2_edge = pydot.Edge(dummy_4, put_together_gap2)
            dot_graph.add_edge(put_together_gap2_edge)
            dummy_5_anchor = pydot.Edge(put_together_gap2, dummy_5)
            dot_graph.add_edge(dummy_5_anchor) 
      
      #fp = open(oppfile, 'a')
      #fp.write('\n======================== Post processing EM stage ========================')
      #fp.close()
      dummy_6 = pydot.Node('5', margin=0, width=0, height=0)
      tree.create_node("I: Post processing EM step", "5", parent="D0")
      dummy_6.set_shape('circle')
      dot_graph.add_node(dummy_6)  
      ret_list = successOrfail_multi(filename, parsing_gaps_5_model_i, "Successfully post processing EM step", "Error in post processing EM step", "Unknown status for post processing EM step")
      count = 0
      for ele in ret_list:
            count = count+1
            #post_process_em = pydot.Node('processing EM step\n%s'%(ele), margin=0, width=0, height=0)
            if  "Successfully" in ele:
                  post_process_em = pydot.Node("I.5.%d"%(count), style='filled', color='green', margin=0, width=0, height=0)
                  #tree.create_node("I.5.%d: %s"%(count,ele), "I.5.%d"%(count), parent="5")
            elif "Error" in ele:
                  post_process_em = pydot.Node("I.5.%d"%(count), style='filled', color='red', margin=0, width=0, height=0)
                  #tree.create_node("I.5.%d: %s"%(count,ele), "I.5.%d"%(count), parent="5")            
            elif "Abscent" in ele:
                  post_process_em = pydot.Node("I.5.%d"%(count), style='filled', color='orange', margin=0, width=0, height=0)
                  #tree.create_node("I.5.%d: %s"%(count,ele), "I.5.%d"%(count), parent="5")
            tree.create_node("I.5.%d: Post processing EM step"%(count), "I.5.%d"%(count), parent="5")
            post_process_em.set_shape('box')
            dot_graph.add_node(post_process_em)
            post_process_em_edge = pydot.Edge(dummy_5, post_process_em)
            dot_graph.add_edge(post_process_em_edge)
            dummy_6_anchor = pydot.Edge(post_process_em, dummy_6)
            dot_graph.add_edge(dummy_6_anchor) 

      #fp = open(oppfile, 'a')
      #fp.write('\n======================== Post processing water refinement preparation ========================')
      #fp.close()
      dummy_7 = pydot.Node('6', margin=0, width=0, height=0)
      tree.create_node("J: File conversion for refinement", "6", parent="D0")
      dummy_7.set_shape('circle')
      dot_graph.add_node(dummy_7)       
      ret_list = successOrfail_multi(filename, parsing_gaps_6_model_i, "Successfully converted for water refinement step", "Error in convertion for water refinement step", "Unknown status for convertion for water refinement step")
      count = 0
      for ele in ret_list:
           count = count+1
           #post_process_wf1 = pydot.Node('converted for water refinement step\n%s'%(ele), margin=0, width=0, height=0)
           if  "Successfully" in ele:
                  post_process_wf1 = pydot.Node("J.6.%d"%(count), style='filled', color='green', margin=0, width=0, height=0)
                  #tree.create_node("J.6.%d: %s"%(count,ele), "J.6.%d"%(count), parent="6")
           elif "Error" in ele:
                  post_process_wf1 = pydot.Node("J.6.%d"%(count), style='filled', color='red', margin=0, width=0, height=0)
                  #tree.create_node("J.6.%d: %s"%(count,ele), "J.6.%d"%(count), parent="6")            
           elif "Abscent" in ele:
                  post_process_wf1 = pydot.Node("J.6.%d"%(count), style='filled', color='orange', margin=0, width=0, height=0)
                  #tree.create_node("J.6.%d: %s"%(count,ele), "J.6.%d"%(count), parent="6")
           tree.create_node("J.6.%d: Conversion for water-refinement"%(count), "J.6.%d"%(count), parent="6")            
           post_process_wf1.set_shape('box')
           dot_graph.add_node(post_process_wf1)
           post_process_wf1_edge = pydot.Edge(dummy_6, post_process_wf1)
           dot_graph.add_edge(post_process_wf1_edge)
           dummy_7_anchor = pydot.Edge(post_process_wf1, dummy_7)
           dot_graph.add_edge(dummy_7_anchor)      

      #fp = open(oppfile, 'a')
      #fp.write('\n======================== Post processing water refinement ========================')
      #fp.close()
      dummy_8 = pydot.Node('7', margin=0, width=0, height=0)
      tree.create_node("k: Post process water refinement", "7", parent="D0")
      dummy_8.set_shape('circle')
      dot_graph.add_node(dummy_8)  
      ret_list = successOrfail_multi(filename, parsing_gaps_7_model_i, "Successfully together running water refinement step in post processing", "Error in running water refinement step in\n post processing", "Unknown status running water refinement step in post processing")
      count = 0
      for ele in ret_list:
           count = count+1
           #post_process_wf2 = pydot.Node('water refinement step\n%s'%(ele), margin=0, width=0, height=0)
           if  "Successfully" in ele:
                  post_process_wf2 = pydot.Node("k.7.%d"%(count), style='filled', color='green', margin=0, width=0, height=0)
                  #tree.create_node("k.7.%d: %s"%(count,ele), "k.7.%d"%(count), parent="7")
           elif "Error" in ele:
                  post_process_wf2 = pydot.Node("k.7.%d"%(count), style='filled', color='red', margin=0, width=0, height=0)
                  #tree.create_node("k.7.%d: %s"%(count,ele), "k.7.%d"%(count), parent="7")            
           elif "Abscent" in ele:
                  post_process_wf2 = pydot.Node("k.7.%d"%(count), style='filled', color='orange', margin=0, width=0, height=0)
                  #tree.create_node("k.7.%d: %s"%(count,ele), "k.7.%d"%(count), parent="7")                       
           tree.create_node("k.7.%d: water-refinment step"%(count), "k.7.%d"%(count), parent="7")
           post_process_wf2.set_shape('box')
           dot_graph.add_node(post_process_wf2)
           post_process_wf2_edge = pydot.Edge(dummy_7, post_process_wf2)
           dot_graph.add_edge(post_process_wf2_edge)
           dummy_8_anchor = pydot.Edge(post_process_wf2, dummy_8)
           dot_graph.add_edge(dummy_8_anchor)  

      dot_graph.write_svg('%s.svg'%(oppfile))
      dot_graph.write_svg('%s.ps2'%(oppfile))   
      print(oppfile)   
      #SVG('%s.svg'%(oppfile))
      dot_graph.write_png('%s.png'%(oppfile))
      tree.to_graphviz('%s.abbr.dot'%(oppfile))
      #tree.show()
      #subprocess.call(["dot", "-Tpng", '%s.abbr.dot'%(oppfile), "-o", '%s.abbr.png'%(oppfile)])
      #tree.to_json('%s.abbr.json'%(oppfile))
      f_abr = open('%s.abbr.txt'%(oppfile), 'w')
      f_abr.write('Green: Success\tRed: Error\tOrange: Undetermined status\n\n')
      f_abr.close()
      tree.save2file('%s.abbr.txt'%(oppfile))
      # txt as pdf
      #pdf = FPDF()
      #pdf.add_page()
      #pdf.add_font('DejaVu', '', 'DejaVuSansCondensed.ttf', uni=True)
      #pdf.set_font('DejaVu', '', 14)      
      #with open('%s.abbr.txt'%(oppfile), 'r') as f:
      #      for x in f:
      #            print(x)
      #            pdf.cell(200, 10, txt = x, ln = 1, align = 'c')

      #pdf.output('%s.abbr.pdf'%(oppfile))



def parseMessageFiles_expert_withparams(filename, oppfile):
      '''
      same as parseMessageFiles_archived
      '''
      print("do-1")

def parseMessageFiles_expert_simple(filename, oppfile):
      '''

      '''
      
      dot_graph = pydot.Dot(graph_type='digraph')
      tree = Tree()
      tree.create_node("DREAM","D0")

      #fp = open(oppfile, 'w')
      #fp.write('\n======================== Prediction stage ========================')
      #fp.close()
      pred_stage = pydot.Node('Prediction')
      pred_stage.set_shape('box3d')
      dot_graph.add_node(pred_stage)
      parsing_for_prediction = "preprediction file"
      ret_list = successOrfail_multi(filename, parsing_for_prediction, "Successfully created file for prediction", "Error in creating file for prediction", "Unknown status for creating file for prediction")
      msgs = "All predicted"
      for ele in ret_list:
            msgs = '%s %s'%(msgs,ele)           
      #pred_stage_i = pydot.Node(msgs, margin=0, width=0, height=0)
      pred_stage_i = pydot.Node('A', margin=0, width=0, height=0)
      tree.create_node("A: Parsing parameter file and predicting %s"%(msgs),"A", parent="D0")      
      pred_stage_i.set_shape('box')
      dot_graph.add_node(pred_stage_i)
      pred_stage_edge = pydot.Edge(pred_stage, pred_stage_i)
      pred_stage_edge.set_label('')
      dot_graph.add_edge(pred_stage_edge)

      #fp = open(oppfile, 'a')
      #fp.write('\n======================== Dendrogram generation stage ========================')
      #fp.close()
      parsing_for_dendrogram = "dendrogram file"
      ret_list = successOrfail_one(filename, parsing_for_dendrogram, "Successfully generated dendrogram", "Error in creating dendrogram", "Unknown status in creating dendrogram")
      msgs = "Dendrograms %s"%(ret_list[0])      
      #dendrogam_stage_i = pydot.Node(msgs, margin=0, width=0, height=0)
      if  "Successfully" in ret_list[0]:
            dendrogam_stage_i = pydot.Node('B', style='filled', color='green', margin=0, width=0, height=0)
            #tree.create_node("B: %s"%(ret_list[0]),"B", parent="D0")
      elif "Error" in ret_list[0]:
            dendrogam_stage_i = pydot.Node('B', style='filled', color='red', margin=0, width=0, height=0)
            #tree.create_node("B: %s"%(ret_list[0]),"B", parent="D0")            
      elif "Abscent" in ret_list[0]:
            dendrogam_stage_i = pydot.Node('B', style='filled', color='orange', margin=0, width=0, height=0)
            #tree.create_node("B: %s"%(ret_list[0]),"B", parent="D0")      
      tree.create_node("B: %s"%(msgs),"B", parent="A")
      dendrogam_stage_i.set_shape('box')
      dot_graph.add_node(dendrogam_stage_i)
      dendrogram_stage_edge = pydot.Edge(pred_stage_i, dendrogam_stage_i)
      dendrogram_stage_edge.set_label('')
      dot_graph.add_edge(dendrogram_stage_edge)      

      #fp = open(oppfile, 'a')
      #fp.write('\n======================== Anchored localize and Gap correct stage ========================')
      #fp.close()
      parsing_for_gap_corr = "gap correct file"
      params_ran = returnMatchingFIlenames(filename, parsing_for_gap_corr)
      ret_list = successOrfail_multi(filename, parsing_for_prediction, "Successfully completed gap correct ", "Error in gap correct module ", "Unknown status for gap correct module ")

      dummy_1 = pydot.Node('+', margin=0, width=0, height=0)
      dummy_1.set_shape('box')
      dot_graph.add_node(dummy_1)

      #for ele in ret_list:      
      #      anchor_and_gapfill = pydot.Node(ele, margin=0, width=0, height=0)
      #      anchor_and_gapfill.set_shape('box')
      #      dot_graph.add_node(anchor_and_gapfill)
      #      anchor_and_gapfill_edge = pydot.Edge(dendrogam_stage_i, anchor_and_gapfill)
      #      anchor_and_gapfill_edge.set_label('')
      #      dot_graph.add_edge(anchor_and_gapfill_edge)

      #      dummy_1_edge = pydot.Edge(anchor_and_gapfill, dummy_1)
      #      dummy_1_edge.set_label('')
      #      dot_graph.add_edge(dummy_1_edge)
      
      dummy_1_edge = pydot.Edge(dendrogam_stage_i, dummy_1)
      dot_graph.add_edge(dummy_1_edge)
      tree.create_node("C:","C", parent="B")

      count_p = 0
      for params in params_ran:
            count_p = count_p+1
            tree.create_node("D%d"%(count_p),"D%d"%(count_p), parent="B")

            #fp = open(oppfile, 'a')
            #fp.write('\n======================== fill gaps for anchored localization ========================')
            #fp.close()
            msgs = "fill gaps and anchored localization %s"%(params)
            #anchor_msg = pydot.Node(msgs, margin=0, width=0, height=0)
            anchor_msg = pydot.Node('E_%d'%(count_p), style='filled', color='orange', margin=0, width=0, height=0)
            if  "Successfully" in ret_list[0]:
                   anchor_msg = pydot.Node('E_%d'%(count_p), style='filled', color='green', margin=0, width=0, height=0)
                   #tree.create_node("B: %s"%(ret_list[0]),"B", parent="D0")
            elif "Error" in ret_list[0]:
                   anchor_msg = pydot.Node('E_%d'%(count_p), style='filled', color='red', margin=0, width=0, height=0)
                   #tree.create_node("B: %s"%(ret_list[0]),"B", parent="D0")            
            elif "Abscent" in ret_list[0]:
                   anchor_msg = pydot.Node('E_%d'%(count_p), style='filled', color='orange', margin=0, width=0, height=0)
                   #tree.create_node("B: %s"%(ret_list[0]),"B", parent="D0")      
            tree.create_node("E_%d: %s"%(count_p,msgs),"E_%d"%(count_p), parent="D%d"%(count_p)) 
            anchor_msg.set_shape('box')
            dot_graph.add_node(anchor_msg)
            anchor_msg_edge = pydot.Edge(dummy_1, anchor_msg)
            anchor_msg_edge.set_label('')
            dot_graph.add_edge(anchor_msg_edge)

            created_filename = "%s"%(params)
            parsing_for_modeller_for_anchored = "initial file for anchored localization %s"%(created_filename)
            ret_list = successOrfail_one(filename, parsing_for_modeller_for_anchored, "Successfully created initial file for anchored localization", "Error in creating initial file for anchored localization", "Unknown status for creating initial file for anchored localization")
            #anchor_ini = pydot.Node('initial\n%s\n%s'%(ret_list[0],params), margin=0, width=0, height=0)
            anchor_ini = pydot.Node('F_%d'%(count_p), style='filled', color='orange', margin=0, width=0, height=0)
            if  "Successfully" in ret_list[0]:
                   anchor_ini = pydot.Node('F_%d'%(count_p), style='filled', color='green', margin=0, width=0, height=0)
                   #tree.create_node("B: %s"%(ret_list[0]),"B", parent="D0")
            elif "Error" in ret_list[0]:
                   anchor_ini = pydot.Node('F_%d'%(count_p), style='filled', color='red', margin=0, width=0, height=0)
                   #tree.create_node("B: %s"%(ret_list[0]),"B", parent="D0")            
            elif "Abscent" in ret_list[0]:
                   anchor_ini = pydot.Node('F_%d'%(count_p), style='filled', color='orange', margin=0, width=0, height=0)
                   #tree.create_node("B: %s"%(ret_list[0]),"B", parent="D0")      
            tree.create_node("F_%d: %s %s"%(count_p,ret_list[0],params),"F_%d"%(count_p), parent="D%d"%(count_p)) 
            anchor_ini.set_shape('box')
            dot_graph.add_node(anchor_ini)
            anchor_ini_edge = pydot.Edge(anchor_msg, anchor_ini)
            anchor_ini_edge.set_label('')
            dot_graph.add_edge(anchor_ini_edge)
         
            #fp = open(oppfile, 'a')
            #fp.write('\n======================== anchored localization stage ========================')
            #fp.close()            
            dummy_2 = pydot.Node('+\n%s'%(params), margin=0, width=0, height=0)
            #tree.create_node('G_%d: anchored localization stage'%(count_p), 'G_%d'%(count_p), parent="D%d"%(count_p))            
            dummy_2.set_shape('box')
            dot_graph.add_node(dummy_2)
            created_filename = "stage2_%s"%(params)
            parsing_for_anchored = "anchored localization file %s"%(created_filename)
            ret_list = successOrfail_multi(filename, parsing_for_anchored, "Successfully completed anchored localization", "Error in anchored localization ", "Unknown status for anchored localization ")
            count_tmp = 0
            for ele in ret_list:
                  count_tmp = count_tmp+1
                  #anchor_loc = pydot.Node('anchored localization\n%s\n%s'%(ele,params), margin=0, width=0, height=0)
                  anchor_loc = pydot.Node('F_%d_%d'%(count_p,count_tmp), style='filled', color='orange', margin=0, width=0, height=0)
                  if  "Successfully" in ret_list[0]:
                        anchor_loc = pydot.Node('F_%d_%d'%(count_p,count_tmp), style='filled', color='green', margin=0, width=0, height=0)
                        #tree.create_node("B: %s"%(ret_list[0]),"B", parent="D0")
                  elif "Error" in ret_list[0]:
                        anchor_loc = pydot.Node('F_%d_%d'%(count_p,count_tmp), style='filled', color='red', margin=0, width=0, height=0)
                        #tree.create_node("B: %s"%(ret_list[0]),"B", parent="D0")            
                  elif "Abscent" in ret_list[0]:
                        anchor_loc = pydot.Node('F_%d_%d'%(count_p,count_tmp), style='filled', color='orange', margin=0, width=0, height=0)
                        #tree.create_node("B: %s"%(ret_list[0]),"B", parent="D0")      
                  tree.create_node("F_%d_%d: %s %s"%(count_p,count_tmp,ele,params),"F_%d_%d"%(count_p,count_tmp), parent="F_%d"%(count_p))                   
                  anchor_loc.set_shape('box')
                  dot_graph.add_node(anchor_loc)
                  anchor_loc_edge = pydot.Edge(anchor_ini, anchor_loc)
                  dot_graph.add_edge(anchor_loc_edge)
                  dummy_2_anchor = pydot.Edge(anchor_loc, dummy_2)
                  dot_graph.add_edge(dummy_2_anchor)                  

            #fp = open(oppfile, 'a')
            #fp.write('\n======================== Ending stage for anchored localization ========================')
            #fp.close()            
            parsing_for_anchored_end = "end of anchored localization %s"%(params)
            ret_list = successOrfail_one(filename, parsing_for_anchored_end, "Successfully ran end stages of anchored localization", "Error in running end stages of anchored localization", "Unknown status for end stages of anchored localization")
            #anchor_end = pydot.Node('end anchored localization\n%s\n%s'%(ret_list[0],params), margin=0, width=0, height=0)
            anchor_end = pydot.Node('G_%d'%(count_p), style='filled', color='orange', margin=0, width=0, height=0)
            if  "Successfully" in ret_list[0]:
                  anchor_end = pydot.Node('G_%d'%(count_p), style='filled', color='green', margin=0, width=0, height=0)
                  #tree.create_node("B: %s"%(ret_list[0]),"B", parent="D0")
            elif "Error" in ret_list[0]:
                  anchor_end = pydot.Node('G_%d'%(count_p), style='filled', color='red', margin=0, width=0, height=0)
                  #tree.create_node("B: %s"%(ret_list[0]),"B", parent="D0")            
            elif "Abscent" in ret_list[0]:
                  anchor_end = pydot.Node('G_%d'%(count_p), style='filled', color='orange', margin=0, width=0, height=0)
                  #tree.create_node("B: %s"%(ret_list[0]),"B", parent="D0")            
            tree.create_node('G_%d: end of anchored localization %s'%(count_p,params), 'G_%d'%(count_p), parent="D%d"%(count_p))
            anchor_end.set_shape('box')
            dot_graph.add_node(anchor_end)
            anchor_end_edge = pydot.Edge(dummy_2, anchor_end)
            dot_graph.add_edge(anchor_end_edge)            

            #fp = open(oppfile, 'a')
            #fp.write('\n======================== Structure generation for gaps ========================')
            #fp.close()   
            dummy_3 = pydot.Node('1\n%s'%(params), margin=0, width=0, height=0)            
            dummy_3.set_shape('box')
            dot_graph.add_node(dummy_3)            
            created_filename = "%s_"%(params)
            parsing_for_gap_start = "structure generation for gaps %s_"%(params)
            tree.create_node('H_%d: structure for gaps'%(count_p),'H_%d'%(count_p), parent="D%d"%(count_p))
            ret_list = successOrfail_multi(filename, parsing_for_gap_start, "Successfully generated starting structure for the gap correct ", "Error in generating starting structure for the gap correct ", "Unknown status in generating starting structure for the gap correct ")
            count_tmp = 0
            for ele in ret_list:
                  count_tmp = count_tmp+1
                  #struct_for_gaps = pydot.Node('structures for gaps\n%s\n%s'%(ele,params), margin=0, width=0, height=0)
                  struct_for_gaps = pydot.Node('H_%d_%d'%(count_p,count_tmp), style='filled', color='orange', margin=0, width=0, height=0)
                  if  "Successfully" in ret_list[0]:
                        struct_for_gaps = pydot.Node('H_%d_%d'%(count_p,count_tmp), style='filled', color='green', margin=0, width=0, height=0)
                        #tree.create_node("B: %s"%(ret_list[0]),"B", parent="D0")
                  elif "Error" in ret_list[0]:
                        struct_for_gaps = pydot.Node('H_%d_%d'%(count_p,count_tmp), style='filled', color='red', margin=0, width=0, height=0)
                        #tree.create_node("B: %s"%(ret_list[0]),"B", parent="D0")            
                  elif "Abscent" in ret_list[0]:
                        struct_for_gaps = pydot.Node('H_%d_%d'%(count_p,count_tmp), style='filled', color='orange', margin=0, width=0, height=0)
                        #tree.create_node("B: %s"%(ret_list[0]),"B", parent="D0")      
                  tree.create_node("H_%d_%d: %s %s"%(count_p,count_tmp,ele,params),"H_%d_%d"%(count_p,count_tmp), parent="H_%d"%(count_p))                   
                  struct_for_gaps.set_shape('box')
                  dot_graph.add_node(struct_for_gaps)
                  struct_for_gaps_edge = pydot.Edge(anchor_end, struct_for_gaps)
                  dot_graph.add_edge(struct_for_gaps_edge)
                  dummy_3_gaps = pydot.Edge(struct_for_gaps, dummy_3)
                  dot_graph.add_edge(dummy_3_gaps)

            #fp = open(oppfile, 'a')
            #fp.write('\n======================== Choosing structure for the gaps ========================')
            #fp.close()   
            dummy_4 = pydot.Node('2\n%s'%(params), margin=0, width=0, height=0)
            dummy_4.set_shape('box')
            dot_graph.add_node(dummy_4) 
            parsing_for_choose_start = "choosing of starting structures for gap regions %s_"%(params)
            tree.create_node('I_%d: choosing of starting structure for gap regions'%(count_p),'I_%d'%(count_p), parent="D%d"%(count_p))
            ret_list = successOrfail_multi(filename, parsing_for_choose_start, "Successfully choosen starting structure for the gap correct ", "Error in choosing starting structure for the gap correct ", "Unknown status in choosing starting structure for the gap correct ")
            count_tmp = 0
            for ele in ret_list:
                  count_tmp = count_tmp+1
                  #choose_struct_for_gaps = pydot.Node('Choose structures for gaps\n%s\n%s'%(ele,params), margin=0, width=0, height=0)
                  choose_struct_for_gaps = pydot.Node('I_%d_%d'%(count_p,count_tmp), style='filled', color='orange', margin=0, width=0, height=0)
                  if  "Successfully" in ret_list[0]:
                        choose_struct_for_gaps = pydot.Node('I_%d_%d'%(count_p,count_tmp), style='filled', color='green', margin=0, width=0, height=0)
                        #tree.create_node("B: %s"%(ret_list[0]),"B", parent="D0")
                  elif "Error" in ret_list[0]:
                        choose_struct_for_gaps = pydot.Node('I_%d_%d'%(count_p,count_tmp), style='filled', color='red', margin=0, width=0, height=0)
                        #tree.create_node("B: %s"%(ret_list[0]),"B", parent="D0")            
                  elif "Abscent" in ret_list[0]:
                        choose_struct_for_gaps = pydot.Node('I_%d_%d'%(count_p,count_tmp), style='filled', color='orange', margin=0, width=0, height=0)
                        #tree.create_node("B: %s"%(ret_list[0]),"B", parent="D0")      
                  tree.create_node("I_%d_%d: %s %s"%(count_p,count_tmp,ele,params),"I_%d_%d"%(count_p,count_tmp), parent="I_%d"%(count_p))                   
                  choose_struct_for_gaps.set_shape('box')
                  dot_graph.add_node(choose_struct_for_gaps)
                  choose_struct_for_gaps_edge = pydot.Edge(dummy_3, choose_struct_for_gaps)
                  dot_graph.add_edge(choose_struct_for_gaps_edge)
                  dummy_4_gaps = pydot.Edge(choose_struct_for_gaps, dummy_4)
                  dot_graph.add_edge(dummy_4_gaps)

            #fp = open(oppfile, 'a')
            #fp.write('\n======================== Putting together chosen structure for the gaps initial ========================')
            #fp.close()   
            dummy_5 = pydot.Node('3\n%s'%(params), margin=0, width=0, height=0)
            dummy_5.set_shape('box')
            dot_graph.add_node(dummy_5)             
            parsing_for_together_gaps = "together choosing starting structures for gap regions %s_"%(params)
            tree.create_node('J_%d: together choosing structure for gap regions'%(count_p),'J_%d'%(count_p), parent="D%d"%(count_p))
            ret_list = successOrfail_multi(filename, parsing_for_together_gaps, "Successfully choosing one off the starting structure for gap correct", "Error in choosing one off the starting structure for gap correct", "Unknown status in choosing one off the starting structure for gap correct")
            count_tmp = 0
            for ele in ret_list:
                  count_tmp = count_tmp+1
                  #put_struct_for_gaps = pydot.Node('Put structures for gaps initial\n%s\n%s'%(ele,params), margin=0, width=0, height=0)
                  put_struct_for_gaps = pydot.Node('J_%d_%d'%(count_p,count_tmp), style='filled', color='orange', margin=0, width=0, height=0)
                  if  "Successfully" in ret_list[0]:
                        put_struct_for_gaps = pydot.Node('J_%d_%d'%(count_p,count_tmp), style='filled', color='green', margin=0, width=0, height=0)
                        #tree.create_node("B: %s"%(ret_list[0]),"B", parent="D0")
                  elif "Error" in ret_list[0]:
                        put_struct_for_gaps = pydot.Node('J_%d_%d'%(count_p,count_tmp), style='filled', color='red', margin=0, width=0, height=0)
                        #tree.create_node("B: %s"%(ret_list[0]),"B", parent="D0")            
                  elif "Abscent" in ret_list[0]:
                        put_struct_for_gaps = pydot.Node('J_%d_%d'%(count_p,count_tmp), style='filled', color='orange', margin=0, width=0, height=0)
                        #tree.create_node("B: %s"%(ret_list[0]),"B", parent="D0")      
                  tree.create_node("J_%d_%d: %s %s"%(count_p,count_tmp,ele,params),"J_%d_%d"%(count_p,count_tmp), parent="J_%d"%(count_p))                   
                  put_struct_for_gaps.set_shape('box')
                  dot_graph.add_node(put_struct_for_gaps)
                  put_struct_for_gaps_edge = pydot.Edge(dummy_4, put_struct_for_gaps)
                  dot_graph.add_edge(put_struct_for_gaps_edge)
                  dummy_5_gaps = pydot.Edge(put_struct_for_gaps, dummy_5)
                  dot_graph.add_edge(dummy_5_gaps)            

            #fp = open(oppfile, 'a')
            #fp.write('\n======================== Putting together chosen structure for the gaps final ========================')
            #fp.close()   
            dummy_6 = pydot.Node('4\n%s'%(params), margin=0, width=0, height=0)
            dummy_6.set_shape('box')
            dot_graph.add_node(dummy_6)             
            parsing_for_gaps_whichchoose = "together chosen structrue for the gaps %s_"%(params)
            tree.create_node('K_%d: put together chosen structure for gaps'%(count_p),'K_%d'%(count_p), parent="D%d"%(count_p))
            ret_list = successOrfail_multi(filename, parsing_for_gaps_whichchoose, "Successfully attached corrected structure for gap correct ", "Error in attaching attached corrected structure for gap correct ", "Unknown status for attaching attached corrected structure for gap correct ")
            count_tmp = 0
            for ele in ret_list:
                  count_tmp = count_tmp+1
                  #put_struct_for_gaps_final = pydot.Node('Put structures for gaps final\n%s\n%s'%(ele,params), margin=0, width=0, height=0)
                  put_struct_for_gaps_final = pydot.Node('K_%d_%d'%(count_p,count_tmp), style='filled', color='orange', margin=0, width=0, height=0)
                  if  "Successfully" in ret_list[0]:
                        put_struct_for_gaps_final = pydot.Node('K_%d_%d'%(count_p,count_tmp), style='filled', color='green', margin=0, width=0, height=0)
                        #tree.create_node("B: %s"%(ret_list[0]),"B", parent="D0")
                  elif "Error" in ret_list[0]:
                        put_struct_for_gaps_final = pydot.Node('K_%d_%d'%(count_p,count_tmp), style='filled', color='red', margin=0, width=0, height=0)
                        #tree.create_node("B: %s"%(ret_list[0]),"B", parent="D0")            
                  elif "Abscent" in ret_list[0]:
                        put_struct_for_gaps_final = pydot.Node('K_%d_%d'%(count_p,count_tmp), style='filled', color='orange', margin=0, width=0, height=0)
                        #tree.create_node("B: %s"%(ret_list[0]),"B", parent="D0")      
                  tree.create_node("K_%d_%d: %s %s"%(count_p,count_tmp,ele,params),"K_%d_%d"%(count_p,count_tmp), parent="K_%d"%(count_p))                   
                  put_struct_for_gaps_final.set_shape('box')
                  dot_graph.add_node(put_struct_for_gaps_final)
                  put_struct_for_gaps_final_edge = pydot.Edge(dummy_5, put_struct_for_gaps_final)
                  dot_graph.add_edge(put_struct_for_gaps_final_edge)
                  dummy_6_gaps = pydot.Edge(put_struct_for_gaps_final, dummy_6)
                  dot_graph.add_edge(dummy_6_gaps)


            #fp = open(oppfile, 'a')
            #fp.write('\n======================== Post processing EM step ========================')
            #fp.close()               
            dummy_7 = pydot.Node('5\n%s'%(params), margin=0, width=0, height=0)
            dummy_7.set_shape('box')
            dot_graph.add_node(dummy_7)             
            parsing_for_em = "post processing EM step %s_"%(params)
            tree.create_node('L_%d: post processing EM step'%(count_p),'L_%d'%(count_p), parent="D%d"%(count_p))
            ret_list = successOrfail_multi(filename, parsing_for_em, "Successfully completed EM steps in post processing", "Error in EM steps in post processing", "Unknown status for EM steps in post processing")     
            count_tmp = 0
            for ele in ret_list:
                  count_tmp = count_tmp+1
                  #post_process_em = pydot.Node('Post processing EM step\n%s\n%s'%(ele,params), margin=0, width=0, height=0)
                  post_process_em = pydot.Node('L_%d_%d'%(count_p,count_tmp), style='filled', color='orange', margin=0, width=0, height=0)
                  if  "Successfully" in ret_list[0]:
                        post_process_em = pydot.Node('L_%d_%d'%(count_p,count_tmp), style='filled', color='green', margin=0, width=0, height=0)
                        #tree.create_node("B: %s"%(ret_list[0]),"B", parent="D0")
                  elif "Error" in ret_list[0]:
                        post_process_em = pydot.Node('L_%d_%d'%(count_p,count_tmp), style='filled', color='red', margin=0, width=0, height=0)
                        #tree.create_node("B: %s"%(ret_list[0]),"B", parent="D0")            
                  elif "Abscent" in ret_list[0]:
                        post_process_em = pydot.Node('L_%d_%d'%(count_p,count_tmp), style='filled', color='orange', margin=0, width=0, height=0)
                        #tree.create_node("B: %s"%(ret_list[0]),"B", parent="D0")      
                  tree.create_node("L_%d_%d: %s %s"%(count_p,count_tmp,ele,params),"L_%d_%d"%(count_p,count_tmp), parent="L_%d"%(count_p)) 
                  post_process_em.set_shape('box')
                  dot_graph.add_node(post_process_em)
                  post_process_em_edge = pydot.Edge(dummy_6, post_process_em)
                  dot_graph.add_edge(post_process_em_edge)
                  dummy_7_gaps = pydot.Edge(post_process_em, dummy_7)
                  dot_graph.add_edge(dummy_7_gaps)

            #fp = open(oppfile, 'a')
            #fp.write('\n======================== Post processing water-refinement step-preparation ========================')
            #fp.close()               
            dummy_8 = pydot.Node('6\n%s'%(params), margin=0, width=0, height=0)
            dummy_8.set_shape('box')
            dot_graph.add_node(dummy_8)             
            parsing_for_convert_waterref = "for water refinement step md_related_gapcorrect_v2/%s_"%(params)
            tree.create_node('M_%d: water refinement conversion'%(count_p),'M_%d'%(count_p), parent="D%d"%(count_p))
            ret_list = successOrfail_multi(filename, parsing_for_convert_waterref, "Successfully converted file for water refinement", "Error in converting file for water refinement", "Unknown status for converting file for water refinement")
            count_tmp = 0
            for ele in ret_list:
                  count_tmp = count_tmp+1
                  #post_process_wr_1 = pydot.Node('Post processing water-refinement\n step-preparation\n%s\n%s'%(ele,params), margin=0, width=0, height=0)
                  post_process_wr_1 = pydot.Node('M_%d_%d'%(count_p,count_tmp), style='filled', color='orange', margin=0, width=0, height=0)
                  if  "Successfully" in ret_list[0]:
                        post_process_wr_1 = pydot.Node('M_%d_%d'%(count_p,count_tmp), style='filled', color='green', margin=0, width=0, height=0)
                        #tree.create_node("B: %s"%(ret_list[0]),"B", parent="D0")
                  elif "Error" in ret_list[0]:
                        post_process_wr_1 = pydot.Node('M_%d_%d'%(count_p,count_tmp), style='filled', color='red', margin=0, width=0, height=0)
                        #tree.create_node("B: %s"%(ret_list[0]),"B", parent="D0")            
                  elif "Abscent" in ret_list[0]:
                        post_process_wr_1 = pydot.Node('M_%d_%d'%(count_p,count_tmp), style='filled', color='orange', margin=0, width=0, height=0)
                        #tree.create_node("B: %s"%(ret_list[0]),"B", parent="D0")      
                  tree.create_node("M_%d_%d: %s %s"%(count_p,count_tmp,ele,params),"M_%d_%d"%(count_p,count_tmp), parent="M_%d"%(count_p))                   
                  post_process_wr_1.set_shape('box')
                  dot_graph.add_node(post_process_wr_1)
                  post_process_wr_1_edge = pydot.Edge(dummy_7, post_process_wr_1)
                  dot_graph.add_edge(post_process_wr_1_edge)
                  dummy_8_gaps = pydot.Edge(post_process_wr_1, dummy_8)
                  dot_graph.add_edge(dummy_8_gaps)            

            #fp = open(oppfile, 'a')
            #fp.write('\n======================== Post processing water-refinement step-preparation ========================')
            #fp.close()               
            dummy_9 = pydot.Node('7\n%s'%(params), margin=0, width=0, height=0)
            dummy_9.set_shape('box')
            dot_graph.add_node(dummy_9)             
            parsing_for_waterref = "water refinement step in post processing %s_"%(params)
            tree.create_node('N_%d: water refinement post processing'%(count_p),'N_%d'%(count_p), parent="D%d"%(count_p))
            ret_list = successOrfail_multi(filename, parsing_for_waterref, "Successfully completed the water refinement step", "Error in the water refinement step", "Unknown status for the water refinement step")
            count_tmp = 0
            for ele in ret_list:
                  count_tmp = count_tmp+1
                  #post_process_wr = pydot.Node('Post processing water-refinement step-preparation %s %s'%(ele,params), margin=0, width=0, height=0)
                  post_process_wr = pydot.Node('N_%d_%d'%(count_p,count_tmp), style='filled', color='orange',margin=0, width=0, height=0)
                  if  "Successfully" in ret_list[0]:
                        post_process_wr = pydot.Node('N_%d_%d'%(count_p,count_tmp), style='filled', color='green', margin=0, width=0, height=0)
                        #tree.create_node("B: %s"%(ret_list[0]),"B", parent="D0")
                  elif "Error" in ret_list[0]:
                        post_process_wr = pydot.Node('N_%d_%d'%(count_p,count_tmp), style='filled', color='red', margin=0, width=0, height=0)
                        #tree.create_node("B: %s"%(ret_list[0]),"B", parent="D0")            
                  elif "Abscent" in ret_list[0]:
                        post_process_wr = pydot.Node('N_%d_%d'%(count_p,count_tmp), style='filled', color='orange', margin=0, width=0, height=0)
                        #tree.create_node("B: %s"%(ret_list[0]),"B", parent="D0")      
                  tree.create_node("N_%d_%d: %s %s"%(count_p,count_tmp,ele,params),"N_%d_%d"%(count_p,count_tmp), parent="N_%d"%(count_p))                   
                  post_process_wr.set_shape('box')
                  dot_graph.add_node(post_process_wr)
                  post_process_wr_edge = pydot.Edge(dummy_8, post_process_wr)
                  dot_graph.add_edge(post_process_wr_edge)
                  dummy_9_gaps = pydot.Edge(post_process_wr, dummy_9)
                  dot_graph.add_edge(dummy_9_gaps)

      dot_graph.write_svg('%s.svg'%(oppfile))
      dot_graph.write_svg('%s.ps2'%(oppfile))
      SVG('%s.svg'%(oppfile))
      dot_graph.write_png('%s.png'%(oppfile))

      #tree.show()
      f_abr = open('%s.abbr.txt'%(oppfile), 'w')
      f_abr.write('Green: Success\tRed: Error\tOrange: Undetermined status\n\n')
      f_abr.close()      
      tree.save2file('%s.abbr.txt'%(oppfile))


def sendMailWithAttachments(recipient, filenames, protname="default", mailbody="attached results"):
       '''
       toaddr: email to which files are being sent
       filename: ['file1','file2',...]
       '''
       port = 587
       server = "smtp-mail.outlook.com"
       sender = "niladrid@iisc.ac.in"
       password = "Niladri#11648"       
       
       # attach body of the mail
       msg = MIMEMultipart()       
       message = "The results attached can also be found in http://pallab.serc.iisc.ernet.in/DREAM/userpdbs/" + protname + "/summary.php"
       msg.attach(MIMEText(message, "plain"))

       #attach the attachments
       for filename in filenames:
             with open(filename, "rb") as fp:
             	lines = fp.readlines()
             	for line in lines:
             		print(line)


def main():
       sender = "niladrid@iisc.ac.in"
       print("To: %s"%(sys.argv[1]))
       print("From: %s"%(sys.argv[2]))
       print("Report file: %s"%(sys.argv[3]))
       print("opp file: %s"%(sys.argv[4]))
       print("flag: %s"%(sys.argv[5]))

       to_mail     = sys.argv[1]
       from_mail   = sys.argv[2]
       report_file = sys.argv[3]
       opp_file    = sys.argv[4]
       flag        = int(sys.argv[5])

       #flag: 1-archive 2-expert simple 3-expert params
       if flag == 1:
          parseMessageFiles_archived(report_file, opp_file)
       elif flag == 2:
          parseMessageFiles_expert_simple(report_file, opp_file)
       elif flag == 3:
          parseMessageFiles_archived(report_file, opp_file)
       

       #parseFilesAndMail(sys.argv[2], sys.argv[1], sys.argv[3])
       
if __name__ == "__main__":
       main()
