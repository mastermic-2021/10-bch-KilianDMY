/* gen is a generator of Fq obtained via ffgen(q) */

[q, t, m] = readvec("input.txt");
\\ q est l'ordre du groupe
\\ le code est t = 35-correcteur

gen = ffgen(q, 'a);

encodefq(i, x) = subst( Pol(digits(i, x.p)), 'x, x );

decodefq(z, x) = fromdigits(Vec(z.pol), x.p);

int2fqx(i, x) = Polrev([ encodefq(z,x) | z<-digits(i, x.p^x.f) ]);

fqx2int(p, x) = fromdigits( [ decodefq(z, x) | z<-Vecrev(p) ], x.p^x.f );



\\ Basile a basé son message sur des erreurs ajoutées dans le code correcteur. C'est donc ce que l'on cherche

syndrome(msg, a, b) = {
  distance = 2*t -1;
  return (sum(i=0, distance, subst(msg, 'x, a^(b+i)) * 'x^i));
}

\\ D'après le cours, avec la méthode de Pade, on a R = S*E avec
\\ R, le polynome message, S le syndrome et E le polynome 'erreur', celui qu'on cherche

\\ Merci Julien pour la fonction 'bestapprPade' qui résout l'équation précedente.

decoder(msg, taille) = {
  my(k, a, b, info, pade, R, E, S_pot, tmp);
	for(k = 0, q,
    a = ffprimroot(gen);
    for(b = 0, q-1,
    	info = List();
      pade = bestapprPade( Mod(syndrome(msg, a, b), x^(2*t)) );
      R = numerator(pade);
      E = denominator(pade);

      for(i = 0, q-2,
      	S_pot = subst(E, 'x, a^(-i));
        if(S_pot == 0, val = subst( (R/deriv(E)) * (x^(b-1)), 'x , a^(-i));
          tmp = fqx2int(val, gen);
          listput(info, tmp)
        )
      );

      if(#info > taille, return (Strchr(Vec(info))));
    );
  );
}

message = int2fqx(m, gen);

secret = decoder(message, 10);
print(secret);
