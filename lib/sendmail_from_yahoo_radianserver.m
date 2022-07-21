function sendmail_from_yahoo_radianserver( receiverAddress , m_subject, m_text )
% SEND an email from radianserver@yahoo.com account. 
% 
% sendmail_from_yahoo_radianserver('l.m.gottwald@amc.nl','Test link','Test message including a <A HREF=<http://www.mathworks.com>>Link</A>');
%


host = 'smtp.mail.yahoo.com';
port  = '465';

setpref( 'Internet','E_mail', 'radianserver@yahoo.com' );
setpref( 'Internet', 'SMTP_Server', host );
setpref( 'Internet', 'SMTP_Username',  'radianserver@yahoo.com'  );
setpref( 'Internet', 'SMTP_Password', 'reconstruction@radian' );

props = java.lang.System.getProperties;
props.setProperty( 'mail.smtp.user', 'radianserver@yahoo.com' );
props.setProperty( 'mail.smtp.host', host );
props.setProperty( 'mail.smtp.port', port );
props.setProperty( 'mail.smtp.starttls.enable', 'true' );
props.setProperty( 'mail.smtp.debug', 'true' );
props.setProperty( 'mail.smtp.auth', 'true' );
props.setProperty( 'mail.smtp.socketFactory.port', port );
props.setProperty( 'mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory' );
props.setProperty( 'mail.smtp.socketFactory.fallback', 'false' );


sendmail( receiverAddress , m_subject, m_text );

end