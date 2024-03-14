  
# libraries to be imported
import sys
from email.mime.base import MIMEBase
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
import smtplib, ssl
from email import encoders

# python front_end/sendMail_v2.py niladri.das.2010@gmail.com /data2/nmr/our_algo_final/protein/2m4k_param_2m4k/sendfiles/sendfiles_list.txt 2m4k

def sendTestMail(recipient):
       '''

       '''
       port = 587
       server = "smtp-mail.outlook.com"
       sender = "dream.nmr@outlook.com"
       password = "dream@nmr"       
       
       # attach body of the mail
       msg = MIMEMultipart()       
       message = "Sending test email"
       msg.attach(MIMEText(message, "plain"))
       SSLcontext = ssl.create_default_context()

       with smtplib.SMTP(server, port) as server:
             server.starttls(context=SSLcontext)
             server.login(sender, password)
             server.sendmail(sender, recipient, msg.as_string())
   
def sendMailWithAttachments(recipient, filenames, protname="default", mailbody="attached results"):
       '''
       toaddr: email to which files are being sent
       filename: ['file1','file2',...]
       '''
       port = 587
       server = "smtp-mail.outlook.com"
       sender = "dream.nmr@outlook.com"
       password = "dream@nmr"       
       
       # attach body of the mail
       msg = MIMEMultipart()       
       message = "The results attached can also be found in http://pallab.serc.iisc.ernet.in/DREAM/userpdbs/" + protname + "/summary.php"
       msg.attach(MIMEText(message, "plain"))

       #attach the attachments
       for filename in filenames:
             with open(filename, "rb") as pdf:
                      attachment = MIMEBase("application", "octet-stream")
                      attachment.set_payload(pdf.read())

             encoders.encode_base64(attachment)

             attachment.add_header(
                      "Content-Disposition",
                     f"attachment; filename= {filename}",
             )
             msg.attach(attachment)

       SSLcontext = ssl.create_default_context()

       with smtplib.SMTP(server, port) as server:
             server.starttls(context=SSLcontext)
             server.login(sender, password)
             server.sendmail(sender, recipient, msg.as_string())
       

def parseFilesAndMail(filenametext, toaddr, protname='default'):
       '''
       
       '''
       with open(filenametext) as f:                
                content = f.readlines()
                
       filenames = [lines.strip() for lines in content]
       #print(filenames)
       sendMailWithAttachments(toaddr, filenames, protname)
       
def main():
       sender = "dream.nmr@outlook.com" # "niladrid@iisc.ac.in"
       print("To: %s"%(sys.argv[1]))
       print("From: %s"%(sender))
       print("attach files: %s"%(sys.argv[2]))
       parseFilesAndMail(sys.argv[2], sys.argv[1], sys.argv[3])
       #sendTestMail(sys.argv[1])

       
if __name__ == "__main__":
       main()       

