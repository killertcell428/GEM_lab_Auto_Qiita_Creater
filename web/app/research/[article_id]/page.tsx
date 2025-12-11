'use client';

import { useEffect, useState } from 'react';
import { useParams, useRouter } from 'next/navigation';
import { api } from '@/lib/api';
import LoadingSpinner from '@/components/LoadingSpinner';
import Link from 'next/link';
import { format } from 'date-fns';
import { ja } from 'date-fns/locale/ja';
import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';

interface ResearchData {
  article_id: string;
  topic?: string;
  research_report?: string;
  plan?: any;
  created_at?: string;
  updated_at?: string;
}

export default function ResearchDetailPage() {
  const params = useParams();
  const router = useRouter();
  const articleId = params.article_id as string;

  const [researchData, setResearchData] = useState<ResearchData | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    if (articleId) {
      loadResearchData();
    }
  }, [articleId]);

  const loadResearchData = async () => {
    try {
      setLoading(true);
      setError(null);
      const data = await api.getArticleResearch(articleId);
      setResearchData(data);
    } catch (err) {
      const errorMessage = err instanceof Error ? err.message : 'リサーチ結果の読み込みに失敗しました';
      setError(errorMessage);
      console.error('リサーチ結果読み込みエラー:', err);
    } finally {
      setLoading(false);
    }
  };

  if (loading) {
    return (
      <div className="flex flex-col items-center justify-center min-h-[400px]">
        <LoadingSpinner />
        <p className="mt-4 text-gray-600 dark:text-gray-400">リサーチ結果を読み込んでいます...</p>
      </div>
    );
  }

  if (error || !researchData) {
    return (
      <div className="p-4">
        <div className="p-4 bg-red-100 border border-red-400 text-red-700 rounded">
          <div className="flex justify-between items-center">
            <span>{error || 'リサーチ結果が見つかりません'}</span>
            <div className="space-x-2">
              <button
                onClick={() => router.push('/research')}
                className="ml-4 px-3 py-1 bg-gray-600 text-white rounded hover:bg-gray-700 text-sm"
              >
                リサーチライブラリに戻る
              </button>
              <button
                onClick={() => router.push('/')}
                className="ml-4 px-3 py-1 bg-blue-600 text-white rounded hover:bg-blue-700 text-sm"
              >
                ダッシュボードに戻る
              </button>
            </div>
          </div>
        </div>
      </div>
    );
  }

  return (
    <div className="space-y-6">
      <div className="flex justify-between items-center">
        <h1 className="text-3xl font-bold text-gray-900 dark:text-white">リサーチ結果</h1>
        <div className="space-x-2">
          <Link
            href={`/articles/${articleId}`}
            className="px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700 transition-colors"
          >
            記事詳細を見る
          </Link>
          <Link
            href="/research"
            className="px-4 py-2 bg-gray-600 text-white rounded-md hover:bg-gray-700 transition-colors"
          >
            リサーチライブラリに戻る
          </Link>
        </div>
      </div>

      <div className="bg-white dark:bg-gray-800 rounded-lg shadow p-6 space-y-6">
        {/* 基本情報 */}
        <div>
          <h2 className="text-xl font-semibold text-gray-900 dark:text-white mb-4">基本情報</h2>
          <div className="space-y-2">
            <div>
              <span className="text-sm font-medium text-gray-600 dark:text-gray-400">記事ID: </span>
              <span className="text-sm text-gray-900 dark:text-white">{researchData.article_id}</span>
            </div>
            {researchData.topic && (
              <div>
                <span className="text-sm font-medium text-gray-600 dark:text-gray-400">トピック: </span>
                <span className="text-sm text-gray-900 dark:text-white">{researchData.topic}</span>
              </div>
            )}
            {researchData.created_at && (
              <div>
                <span className="text-sm font-medium text-gray-600 dark:text-gray-400">作成日: </span>
                <span className="text-sm text-gray-900 dark:text-white">
                  {format(new Date(researchData.created_at), 'yyyy年MM月dd日 HH:mm', { locale: ja })}
                </span>
              </div>
            )}
          </div>
        </div>

        {/* リサーチ結果 */}
        {researchData.research_report && (
          <div>
            <h2 className="text-xl font-semibold text-gray-900 dark:text-white mb-4">リサーチ結果</h2>
            <div className="p-4 bg-gray-50 dark:bg-gray-900 rounded">
              <div className="prose dark:prose-invert max-w-none">
                <ReactMarkdown remarkPlugins={[remarkGfm]}>
                  {researchData.research_report}
                </ReactMarkdown>
              </div>
            </div>
          </div>
        )}

        {/* 記事構成プラン */}
        {researchData.plan && (
          <div>
            <h2 className="text-xl font-semibold text-gray-900 dark:text-white mb-4">記事構成プラン</h2>
            <div className="p-4 bg-gray-50 dark:bg-gray-900 rounded">
              <pre className="text-sm overflow-auto max-h-96">
                {JSON.stringify(researchData.plan, null, 2)}
              </pre>
            </div>
          </div>
        )}

        {!researchData.research_report && !researchData.plan && (
          <div className="text-center py-12">
            <p className="text-gray-600 dark:text-gray-400">リサーチ結果がありません</p>
          </div>
        )}
      </div>
    </div>
  );
}

