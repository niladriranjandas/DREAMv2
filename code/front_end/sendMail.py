  
# libraries to be imported
import sys
import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.base import MIMEBase
from email import encoders

FROM_MAIL = "dream.email.sender@gmail.com"
FROM_PASS = "niladri@123"
   
def sendMailWithAttachments(toaddr, filenames, protname="default", mailbody="attached results"):
       '''
       toaddr: email to which files are being sent
       filename: ['file1','file2',...]
       '''
       # instance of MIMEMultipart
       msg = MIMEMultipart()
  
       # storing the senders email address  
       msg['From'] = FROM_MAIL
  
       # storing the receivers email address 
       msg['To'] = toaddr
  
       # storing the subject 
       msg['Subject'] = "DREAM : " + protname 
  
       # string to store the body of the mail
       #body = "Body_of_the_mail"
       body = "The results attached can also be found in http://pallab.serc.iisc.ernet.in/DREAM/userpdbs/" + protname + "/summary.php"
  
       # attach the body with the msg instance
       msg.attach(MIMEText(body, 'plain'))
       
       # ------- attachment ---------------# 
       # open the file to be sent 
       for filename in filenames:
                print(filename)
                #filename =  "/home/niladri/Downloads/rough/check_recent/2m4k/gaps_no_S/2m4k_param_2m4k_1_ancloc_refine_4_em_nosol_wref.pdb"  # "File_name_with_extension"
                attachment = open(filename, "rb")
  
                # instance of MIMEBase and named as p
                p = MIMEBase('application', 'octet-stream')
  
                # To change the payload into encoded form
                p.set_payload((attachment).read())
  
                # encode into base64
                encoders.encode_base64(p)
   
                p.add_header('Content-Disposition', "attachment; filename= %s" % filename)
  
                # attach the instance 'p' to instance 'msg'
                msg.attach(p)


       # creates SMTP session
       s = smtplib.SMTP('smtp.gmail.com', 587)
  
       # start TLS for security
       s.starttls()
  
       # Authentication
       s.login(FROM_MAIL, FROM_PASS)
  
       # Converts the Multipart msg into a string
       text = msg.as_string()
  
       # sending the mail
       s.sendmail(FROM_MAIL, toaddr, text)
  
       # terminating the session
       s.quit()

def parseFilesAndMail(filenametext, toaddr, protname='default'):
       '''
       
       '''
       with open(filenametext) as f:                
                content = f.readlines()
                
       filenames = [lines.strip() for lines in content]
       #print(filenames)
       sendMailWithAttachments(toaddr, filenames, protname)
       
def main():
       print("To: %s"%(sys.argv[1]))
       print("From: %s"%(FROM_MAIL))
       print("attach files: %s"%(sys.argv[2]))
       parseFilesAndMail(sys.argv[2], sys.argv[1], sys.argv[3])
       
if __name__ == "__main__":
       main()       

