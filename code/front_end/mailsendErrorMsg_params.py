import re
import sys
from email.mime.base import MIMEBase
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
import smtplib, ssl
from email import encoders

# python front_end/sendMail_v2.py niladri.das.2010@gmail.com /data2/nmr/our_algo_final/protein/2m4k_param_2m4k/sendfiles/sendfiles_list.txt 2m4k

def successOrfail_one(filename, oppfile, srch_pattern, success_msg, fail_msg, abscent_msg):
       '''

       '''
       regex = re.compile(srch_pattern)
       with open(filename, 'r') as fp, open(oppfile,'a') as op:
            op.write("\n---%s---"%(srch_pattern))
            for line in fp:
                  k = regex.findall(line)
                  position = False
                  if str(k) != '[]':
                        position = True
                        if "success" in line.lower():
                               #print(line)
                               op.write('\n%s'%(success_msg))
                        elif "error" in line.lower():
                               #print(line)
                               op.write('\n%s'%(fail_msg))
                        else:
                               #print(line)
                               op.write('\n%s'%(abscent_msg))

       #if str(k) == '[]':
       #     op = open(oppfile, 'a')
       #     op.write('\n%s'%(abscent_msg))
       #     op.close()

def successOrfail_multi(filename, oppfile, srch_pattern, success_msg, fail_msg, abscent_msg):
       '''

       '''
       regex = re.compile(srch_pattern)
       c=0
       with open(filename, 'r') as fp, open(oppfile,'a') as op:
            op.write("\n---%s---"%(srch_pattern))
            for line in fp:
                  k = regex.findall(line)
                  position = False
                  if str(k) != '[]':
                        position = True
                        c = c+1
                        if "success" in line.lower():
                               op.write('\n%s for model-%d'%(success_msg,c))
                        elif "error" in line.lower():
                               op.write('\n%s for model-%d'%(fail_msg,c))
                        else:
                               op.write('\n%s for model-%d'%(abscent_msg,c))

       #if str(k) == '[]':
       #     op = open(oppfile, 'a')
       #     op.write('\n%s'%(abscent_msg))
       #     op.close()


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
      
      parsing_file_pat    = "parsing the input paramter files"
      parsing_div_frags   = "dividing and individual modeling section"
      parsing_consolidate = "consolidating the divided parts"
      parsing_gaps_model  = "modelled the gaps post processing steps"
      parsing_gaps_post   = "modelled the gaps post processing steps"

      opfilename = '%s_sendmail.txt'%(filename)

      fp = open(oppfile, 'w')
      fp.write('\n======================== Parsing parameter file ========================')
      fp.close()
      successOrfail_one(filename, oppfile, parsing_file_pat, "Successfully parsed the input file.", "Error while dividing into fragments and their modelling.", "Abscent status while dividing into fragments and their modelling.")
      fp = open(oppfile, 'a')
      fp.write('\n======================== Dividing into fragments and modelling ========================')
      fp.close()
      successOrfail_one(filename, oppfile, parsing_div_frags, "Successfully divided into fragments and modelled them.", "Error while dividing into fragments and their modelling.", "Abscent status while dividing into fragments and their modelling.")
      fp = open(oppfile, 'a')
      fp.write('\n======================== Consolidating structures for fragments ========================')
      fp.close()
      successOrfail_one(filename, oppfile, parsing_consolidate, "Successfully while consolidating structure for fragments.", "Error while consolidating structure for fragments.", "Abscent status while consolidating structure for fragments.")
      fp = open(oppfile, 'a')
      fp.write('\n======================== Modelling the gaps ========================')
      fp.close()
      successOrfail_one(filename, oppfile, parsing_gaps_model, "Successfully modelled tha gaps.", "Error while modelling gaps.", "Abscent status Error while modelling gaps.")
      fp = open(oppfile, 'a')
      fp.write('\n======================== Prediction stage ========================')
      fp.close()
      successOrfail_one(filename, oppfile, parsing_gaps_post, "Successfully modelled tha gaps.", "Error while modelling gaps.", "Abscent status Error while modelling gaps.")      

      parsing_gaps_1_model_i = "generating starting structure generation for gaps"
      parsing_gaps_2_model_i = "in choosing starting structures for gap regions"
      parsing_gaps_3_model_i = "together choosing starting structures for gap regions"
      parsing_gaps_4_model_i = "together chosen structrue for the gaps"
      parsing_gaps_5_model_i = "post processing EM step"
      parsing_gaps_6_model_i = "for water refinement step"
      parsing_gaps_7_model_i = "water refinement step in post processing"

      fp = open(oppfile, 'a')
      fp.write('\n======================== Gap correction starting sturcture ========================')
      fp.close()
      successOrfail_multi(filename, oppfile, parsing_gaps_1_model_i, "Successfully modelled starting structure for gaps", "Error in modelling starting structure for gaps", "Unknown status for modelling starting structure for gaps")
      fp = open(oppfile, 'a')
      fp.write('\n======================== Generating structures for the gap correct ========================')
      fp.close()
      successOrfail_multi(filename, oppfile, parsing_gaps_2_model_i, "Successfully choosen starting structures for gap regions", "Error in choosing starting structures for gap regions", "Unknown status in choosing starting structures for gap regions")
      fp = open(oppfile, 'a')
      fp.write('\n======================== Putting together the gaps corrected initial ========================')
      fp.close()
      successOrfail_multi(filename, oppfile, parsing_gaps_3_model_i, "Successfully put together choosing starting structures for gap regions", "Error in putting together choosing starting structures for gap regions", "Unknown status for putting together choosing starting structures for gap regions")
      fp = open(oppfile, 'a')
      fp.write('\n======================== Putting together the gaps corrected final ========================')
      fp.close()
      successOrfail_multi(filename, oppfile, parsing_gaps_4_model_i, "Successfully together chosen structrue for the gaps", "Error in together chosen structrue for the gaps", "Unknown status for together chosen structrue for the gaps")
      fp = open(oppfile, 'a')
      fp.write('\n======================== Post processing EM stage ========================')
      fp.close()
      successOrfail_multi(filename, oppfile, parsing_gaps_5_model_i, "Successfully post processing EM step", "Error in post processing EM step", "Unknown status for post processing EM step")
      fp = open(oppfile, 'a')
      fp.write('\n======================== Post processing water refinement preparation ========================')
      fp.close()
      successOrfail_multi(filename, oppfile, parsing_gaps_6_model_i, "Successfully converted for water refinement step", "Error in convertion for water refinement step", "Unknown status for convertion for water refinement step")
      fp = open(oppfile, 'a')
      fp.write('\n======================== Post processing water refinement ========================')
      fp.close()
      successOrfail_multi(filename, oppfile, parsing_gaps_7_model_i, "Successfully together running water refinement step in post processing", "Error in running water refinement step in post processing", "Unknown status running water refinement step in post processing")




def parseMessageFiles_expert_withparams(filename, oppfile):
      '''
      same as parseMessageFiles_archived
      '''
      print("do-1")

def parseMessageFiles_expert_simple(filename, oppfile):
      '''

      '''
      

      fp = open(oppfile, 'w')
      fp.write('\n======================== Prediction stage ========================')
      fp.close()
      parsing_for_prediction = "preprediction file"
      successOrfail_multi(filename, oppfile, parsing_for_prediction, "Successfully created file for prediction", "Error in creating file for prediction", "Unknown status for creating file for prediction")

      fp = open(oppfile, 'a')
      fp.write('\n======================== Dendrogram generation stage ========================')
      fp.close()
      parsing_for_dendrogram = "dendrogram file"
      successOrfail_one(filename, oppfile, parsing_for_dendrogram, "Successfully generated dendrogram", "Error in creating dendrogram", " Unknown status in creating dendrogram")

      fp = open(oppfile, 'a')
      fp.write('\n======================== Anchored localize and Gap correct stage ========================')
      fp.close()
      parsing_for_gap_corr = "gap correct file"
      params_ran = returnMatchingFIlenames(filename, parsing_for_gap_corr)
      successOrfail_multi(filename, oppfile, parsing_for_prediction, "Successfully completed gap correct", "Error in gap correct module", "Unknown status for gap correct module")

      print("here")
      for params in params_ran:

            fp = open(oppfile, 'a')
            fp.write('\n======================== fill gaps for anchored localization ========================')
            fp.close()
            created_filename = "%s"%(params)
            parsing_for_modeller_for_anchored = "initial file for anchored localization %s"%(created_filename)
            successOrfail_one(filename, oppfile, parsing_for_modeller_for_anchored, "Successfully created initial file for anchored localization", "Error in creating initial file for anchored localization", "Unknown status for creating initial file for anchored localization")

            fp = open(oppfile, 'a')
            fp.write('\n======================== anchored localization stage ========================')
            fp.close()            
            created_filename = "stage2_%s"%(params)
            parsing_for_anchored = "anchored localization file %s"%(created_filename)
            successOrfail_multi(filename, oppfile, parsing_for_anchored, "Successfully completed anchored localization for ", "Error in anchored localization", "Unknown status for anchored localization")

            fp = open(oppfile, 'a')
            fp.write('\n======================== Ending stage for anchored localization ========================')
            fp.close()
            parsing_for_anchored_end = "end of anchored localization"
            successOrfail_one(filename, oppfile, parsing_for_anchored_end, "Successfully ran end stages of anchored localization", "Error in running end stages of anchored localization", "Unknown status for end stages of anchored localization")
            

            fp = open(oppfile, 'a')
            fp.write('\n======================== Structure generation for gaps ========================')
            fp.close()   
            created_filename = "%s_"%(params)
            parsing_for_gap_start = "structure generation for gaps"
            successOrfail_multi(filename, oppfile, parsing_for_gap_start, "Successfully generated starting structure for the gap correct", "Error in generating starting structure for the gap correct", "Unknown status in generating starting structure for the gap correct")

            fp = open(oppfile, 'a')
            fp.write('\n======================== Choosing structure for the gaps ========================')
            fp.close()   
            parsing_for_choose_start = "choosing starting structures for gap regions"
            successOrfail_multi(filename, oppfile, parsing_for_choose_start, "Successfully choosen starting structure for the gap correct", "Error in choosing starting structure for the gap correct", "Unknown status in choosing starting structure for the gap correct")

            fp = open(oppfile, 'a')
            fp.write('\n======================== Putting together chosen structure for the gaps initial ========================')
            fp.close()   
            parsing_for_together_gaps = "together choosing starting structures for gap regions"
            successOrfail_multi(filename, oppfile, parsing_for_together_gaps, "Successfully choosing one off the starting structure for gap correct", "Error in choosing one off the starting structure for gap correct", "Unknown status in choosing one off the starting structure for gap correct")

            fp = open(oppfile, 'a')
            fp.write('\n======================== Putting together chosen structure for the gaps final ========================')
            fp.close()   
            parsing_for_gaps_whichchoose = "together chosen structrue for the gaps"
            successOrfail_multi(filename, oppfile, parsing_for_gaps_whichchoose, "Successfully attached corrected structure for gap correct", "Error in attaching attached corrected structure for gap correct", "Unknown status for attaching attached corrected structure for gap correct")

            fp = open(oppfile, 'a')
            fp.write('\n======================== Post processing EM step ========================')
            fp.close()               
            parsing_for_em = "post processing EM step"
            successOrfail_multi(filename, oppfile, parsing_for_em, "Successfully completed EM steps in post processing", "Error in EM steps in post processing", "Unknown status for EM steps in post processing")     

            fp = open(oppfile, 'a')
            fp.write('\n======================== Post processing water-refinement step-preparation ========================')
            fp.close()               
            parsing_for_convert_waterref = "for water refinement step"
            successOrfail_multi(filename, oppfile, parsing_for_convert_waterref, "Successfully converted file for water refinement", "Error in converting file for water refinement", "Unknown status for converting file for water refinement")

            fp = open(oppfile, 'a')
            fp.write('\n======================== Post processing water-refinement step-preparation ========================')
            fp.close()               
            parsing_for_waterref = "water refinement step in post processing"
            successOrfail_multi(filename, oppfile, parsing_for_waterref, "Successfully completed the water refinement step", "Error in the water refinement step", "Unknown status for the water refinement step")

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