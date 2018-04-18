function  send_mail_upon_finished(subject, content, recipient)

mailAddress = '18850544602@163.com'; 
password = 'wy18850544602';  
setpref('Internet','E_mail',mailAddress); 
setpref('Internet','SMTP_Server','smtp.163.com');
setpref('Internet','SMTP_Username',mailAddress);
setpref('Internet','SMTP_Password',password); 
props = java.lang.System.getProperties; 
props.setProperty('mail.smtp.auth','true');
sendmail(recipient,subject,content)


end

