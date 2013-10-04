// SVN: $Id: MakeHtmlDoc.C 46 2012-12-07 02:21:19Z mushkar $
{
  gROOT->Reset();
  gSystem->Load( "libFit.so" );
  THtml h;
  h.MakeAll();
  return;
}
