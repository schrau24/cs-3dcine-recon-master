function sendmail_from_textfile(mail_recipient, mail_subject, mail_text_file)
% sendmail_from_amc(mail_recipient, mail_subject, mail_text)
%
% mail_recipient:   email adress            (char)
% mail_subject:     email subject           (char)
% mail_text_file:   location of email body  (char)
%
% Example: 
% sendmail_from_textfile('l.m.gottwald@amc.uva.nl', 'Testmail', '/home/testfile.txt')
%
% 13-Aug-2018 l.m.gottwald@amc.uva.nl & e.p.grootes@amc.uva.nl

error_id = 0;
if ~ischar(mail_recipient); error_id = 1; end
if ~ischar(mail_subject);   error_id = 1; end
if ~ischar(mail_text_file); error_id = 1; end
if error_id; error('Check class of input variables!'); end

system(sprintf('mail -s "%s" %s < "%s"',mail_subject,mail_recipient,mail_text_file));
end