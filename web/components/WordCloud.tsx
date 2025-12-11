'use client';

import { useEffect, useRef } from 'react';

interface WordCloudProps {
  words: Array<{ text: string; weight: number }>;
  width?: number;
  height?: number;
}

export default function WordCloudComponent({ words, width = 800, height = 400 }: WordCloudProps) {
  const canvasRef = useRef<HTMLCanvasElement>(null);

  useEffect(() => {
    if (words.length === 0 || !canvasRef.current) return;

    // wordcloud2の初期化（動的インポート）
    import('wordcloud').then((WordCloudLib) => {
      if (canvasRef.current) {
        WordCloudLib.default(canvasRef.current, {
          list: words,
          gridSize: Math.round(16 * (width / 1024)),
          weightFactor: function (size: number) {
            return Math.pow(size, 2.3) * (width / 1024);
          },
          fontFamily: 'Arial, sans-serif',
          color: function () {
            return '#' + Math.floor(Math.random() * 16777215).toString(16);
          },
          rotateRatio: 0.5,
          rotationSteps: 2,
          backgroundColor: 'transparent'
        });
      }
    }).catch((err) => {
      console.error('WordCloud読み込みエラー:', err);
    });
  }, [words, width, height]);

  if (words.length === 0) {
    return (
      <div className="bg-white dark:bg-gray-800 rounded-lg shadow p-6">
        <div className="flex items-center justify-center h-64">
          <p className="text-gray-600 dark:text-gray-400">キーワードデータがありません</p>
        </div>
      </div>
    );
  }

  return (
    <div className="bg-white dark:bg-gray-800 rounded-lg shadow p-6">
      <h3 className="text-lg font-semibold mb-4">キーワード分析</h3>
      <canvas
        ref={canvasRef}
        width={width}
        height={height}
        className="w-full"
      />
    </div>
  );
}

