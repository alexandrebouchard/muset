package muset;


import briefj.Indexer;



public class Alphabet
{
  public final Indexer<Letter> indexer = new Indexer<Letter>();
  private int letterStringLength = 0;
  
  public int getMaxLetterStringLength() 
  { 
    return letterStringLength; 
  }
  
  public Alphabet() {}
  
  public Alphabet(Alphabet alphabet)
  {
    for (Letter letter : alphabet.indexer.objectsList())
      this.getLetter(letter.contents);
  }
  
  public boolean containsLetter(String letterContents)
  {
    Letter temp = new Letter(letterContents);
    return indexer.containsObject(temp);
  }

  public Letter getLetter(String contents)
  {
    Letter result = new Letter(contents);
    if (!indexer.containsObject(result))
    {
      indexer.addToIndex(result);
      letterStringLength = Math.max(letterStringLength, contents.length());
    }
    return result;
  }
  
  public String toPaddedString(Letter letter)
  {
    return toPaddedString(letter, letterStringLength);
  }
  
  public static String toPaddedString(Letter letter, int letterStringLength)
  {
    String result = letter.contents;
    while (result.length() < letterStringLength)
      result = result + " ";
    return result;
  }
  
  public static class Letter
  {
    private final String contents;

    private Letter(String contents)
    {
      this.contents = contents;
    }
    
    @Override
    public String toString()
    {
      return contents;
    }

    @Override
    public int hashCode()
    {
      final int prime = 31;
      int result = 1;
      result = prime * result + ((contents == null) ? 0 : contents.hashCode());
      return result;
    }

    @Override
    public boolean equals(Object obj)
    {
      if (this == obj)
        return true;
      if (obj == null)
        return false;
      if (getClass() != obj.getClass())
        return false;
      Letter other = (Letter) obj;
      if (contents == null)
      {
        if (other.contents != null)
          return false;
      } else if (!contents.equals(other.contents))
        return false;
      return true;
    }
  }
}
