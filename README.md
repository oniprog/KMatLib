KMatLib
=======

Tiny matrix calculation library.

The features are as follows.

1. Use stack instated of heap. 
2. C++ template library. And simple. (in other words, too limited function and under development ).
3. License is public domain ( my writing part ). (Notice: Obviously Google C++ Unit Testing framework part is another licesence. )

By using stack, high performace can gain compared with ordinary matrix library.
However, stack is so limited resource compared with heap. Sometime stackover flows may be occur. Therefore this library can adopt for tiny matrix( N <= 10 ).
 
日本語での説明
=======
このライブラリは小さな行列に対する計算ライブラリです．

次の特徴を持ちます

1. スタックを用います．(ヒープではなく）
2. C++テンプレートライブラリでシンプルです．（言い換えれば，機能が限定されていて，開発途上です)
3. Public Domainライセンスです(もちろん私が書いた部分だけ)．利便性のためにGoogle C++ Testing Frameworkのソースを含めていますが，それは当然別ライセンスなのでご注意ください)

スタックを使うことで，通常の行列ライブラリよりは高い性能が得られます．しかし，スタックはヒープと比べるとかなり限られた資源です．気をつけないと，スタックオーバーフローが起こる可能性があります．このライブラリは，N<=10くらいの小さな行列に対してのみ適用することができるでしょう．

oniprog  tkuman@gmail.com

